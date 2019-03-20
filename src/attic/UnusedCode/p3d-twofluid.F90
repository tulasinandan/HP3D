!******************************************************************************
!              3D particle code: main module, two-fluid version
!                           Andreas Zeiler, 2001
!
!                       LATEST CHANGE: July 20, 2001
!
!******************************************************************************

#include "param"
#ifndef twofluid
  switch twofluid must be set for twofluid version
#endif

#ifndef D_b
#define D_b 0.
#endif
#ifndef D_pe
#define D_pe 0.
#endif
#ifndef D_pi
#define D_pi 0.
#endif
#ifndef D_j
#define D_j 0.
#endif
#ifndef D_n
#define D_n 0.
#endif
#ifndef fixedsteps
#define fixedsteps 1
#endif
#ifndef GS_LEX
#ifndef GS_RB
#ifndef JACOBI
#define JACOBI
#endif
#endif
#endif

! assign unique number to each set of boundary conditions
#define periodic 0

Program p3d_twofluid
  Use pe_env
  Use partfield
  Use timer
  Implicit none

  Integer :: n,nsteps,i,filecount,x,y,z,pe
  Character*16 :: buf
  Real(kind=drk) :: t_diag
  Real(kind=drk) :: stoptime, presenttime
  Logical :: terminate=.false.
  Integer :: mpi_err
#ifdef sn6311
  Logical, External :: mpp_terminate
#endif
#define movie_variables
#include movie_header
#undef movie_variables

  call init_pe_env()
  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** p3d-twofluid: compiled for different number of PEs *****'
    call exitallpes()
  endif

#ifdef sn6311
! check, if operator requested termination
  if (mpp_terminate()) call MPI_Finalize(mpi_err)
#endif

! determine stop time according to start time and maximum allowed run time
  stoptime=MPI_Wtime()+maxruntime*60
  call MPI_Bcast(stoptime,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)

! read number of start-up file
  call our_getarg(1,buf)
  if (buf == ' ') then
     filecount = 0
  else
     read (buf,*) filecount
  endif
! read number of time steps
  call our_getarg(2,buf)
  if (buf == ' ') then
     nsteps = 1
  else
     read (buf,*) nsteps
  endif

! write parameters to log-file
  if (myproc==0) then
     write(6,*) '********* p3d-twofluid *********'
     write(6,*) '********** parameters **********'
     write(6,*) '     grid points nx = ',nx
     write(6,*) '     grid points ny = ',ny
     write(6,*) '     grid points nz = ',nz
     write(6,*) '     processors pex = ',pex
     write(6,*) '     processors pey = ',pey
     write(6,*) '     processors pez = ',pez
     write(6,*) '     edge length lx = ',lx
     write(6,*) '     edge length ly = ',ly
     write(6,*) '     edge length lz = ',lz
     write(6,*) '     d_e2 = ',d_e2
     write(6,*) '     dt = ',dt
     write(6,*) '     gamma_i = ',gam_i
     write(6,*) '     gamma_e = ',gam_e
#ifdef n_base
     write(6,*) '     n_base = ',n_base
#endif
#ifdef isothermal
     write(6,*) '     isothermal simulation'
     write(6,*) '        T_e = ',T_e
     write(6,*) '        T_i = ',T_i
#endif
#ifdef double_byte
     write(6,*) '     two-byte movie output'
#endif
     write(6,*) '  boundary conditions:'
#if (boundary_condition == periodic)
     write(6,*) '       periodic'
#else
  we never should get here
#endif
     write(6,*) '  fourth-order viscosities:'
     write(6,*) '       D_b  = ',D_b
     write(6,*) '       D_j  = ',D_j
     write(6,*) '       D_n  = ',D_n
     write(6,*) '       D_pe = ',D_pe
     write(6,*) '       D_pi = ',D_pi
     write(6,*) '  multigrid parameters:'
     write(6,*) '       nu1 = ',nu1
     write(6,*) '       nu2 = ',nu2
     write(6,*) '       nu3 = ',nu3
     write(6,*) '       eps = ',eps
#ifdef alt_restrict_prolong
     write(6,*) '       alt_restrict_prolong'
#endif
#ifdef miniter
     write(6,*) '       miniter =',miniter
#endif
#ifdef fixiter
     write(6,*) '       fixiter =',fixiter
#else
     write(6,*) '       fixiter =',6
#endif
#ifdef JACOBI
     write(6,*) '       smoother: Jacobi'
#endif
#ifdef GS_LEX
     write(6,*) '       smoother: GS_LEX'
#endif
#ifdef GS_RB
     write(6,*) '       smoother: GS_RB'
#endif
#ifndef WCYCLE
     write(6,*) '       V-CYCLE'
#else
     write(6,*) '       W-CYCLE'
#endif
     write(6,*) '       fixedsteps         = ', fixedsteps
     write(6,*) '  default precision       = ',drk
     write(6,*) '  precision for particles = ',prk
     write(6,*) '********** parameters **********'
  endif

! read start-up file
  call input(filecount,.true.)
! fix divergence of b
  call divergence(b1x,b1y,b1z)
! initialize b2, j2, pe2, pi2, n2
  do z=1,nz; do y=1,ny; do x=1,nx
    b2x(x,y,z)=b1x(x,y,z); b2y(x,y,z)=b1y(x,y,z); b2z(x,y,z)=b1z(x,y,z)
    j2x(x,y,z)=j1x(x,y,z); j2y(x,y,z)=j1y(x,y,z); j2z(x,y,z)=j1z(x,y,z)
    pe2(x,y,z)=pe1(x,y,z); pi2(x,y,z)=pi1(x,y,z); n2(x,y,z)=n1(x,y,z)
  enddo; enddo; enddo
! intialize movie files, write first frames
#define movie_open
#include movie_header
#undef movie_open
  if (t .eq. 0.) then
#define movie_write
#include movie_header
#undef movie_write
  endif
! main loop
  call init_timer()
  do n=1,nsteps
    call energy()
    call steptwofluid(mod(n,fixedsteps)>0)
    t=t+dt
#define movie_write
#include movie_header
#undef movie_write
! check, if time limit reached
    presenttime=MPI_Wtime()
    call MPI_Bcast(presenttime,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
    if (presenttime .ge. stoptime) then
      if (myproc == 0) write (6,*) 'time limit reached, save data'
      terminate=.true.
    endif
#ifdef sn6311
! check, if operator requested termination
    if (mpp_terminate()) then
      if (myproc == 0) write (6,*) 'terminating on external request'
      terminate=.true.
    endif
#endif
    if (terminate) exit
  enddo

! write output file
  filecount=filecount+1
  call output(filecount,.true.)
  call read_timer()
! close movie-files
  if (myproc == 0) then
#define movie_close
#include movie_header
#undef movie_close
  endif
  call MPI_Finalize(mpi_err)
End Program p3d_twofluid

!-------------------------------------------------------------------------------

Subroutine our_getarg(num,buf)

  Implicit none
  Integer num
  Character*(*) buf
#ifdef CRAY
  Integer ilen, ierr
#endif

  buf = ' '

#ifdef CRAY
  call pxfgetarg(num,buf,ilen,ierr)
#endif
#ifdef AIX
  call getarg(num,buf)
#endif
#ifdef HITACHI
  call getarg(num+1,buf)
#endif

End Subroutine our_getarg
