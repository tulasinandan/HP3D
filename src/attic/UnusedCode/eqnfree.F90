!******************************************************************************
!                      3D particle code: main module
!                           Andreas Zeiler, 1999
!
!                        LATEST CHANGE: July 4, 2001
!
!******************************************************************************

#include "param"
#ifndef skip_poisson
#define skip_poisson 1
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
#define gem_reconnection 1
#define gem_3d 2

Program p3d
  Use pe_env
  Use partfield
  Use timer
  Implicit none

!  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: temp_phi
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
    write(6,*) '***** p3d: compiled for different number of PEs *****'
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
!  call our_getarg(1,buf)
!  if (buf == ' ') then
!     filecount = 0
!  else
!     read (buf,*) filecount
!  endif
! read number of time steps
!  call our_getarg(2,buf)
!  if (buf == ' ') then
!     nsteps = 1
!  else
!     read (buf,*) nsteps
!  endif

   nsteps = NSTEPSMICRO
   filecount = FILECOUNT

! write parameters to log-file
  if (myproc==0) then
     write(6,*) '**********    p3d     **********'
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
     write(6,*) '     c_2 = ',c_2
     write(6,*) '     m_e = ',m_e
     write(6,*) '     dt = ',dt
     write(6,*) '     substeps = ',substeps
#ifdef fourth_order
     write(6,*) '     fourth order operators'
#endif
#ifdef double_byte
     write(6,*) '     two-byte movie output'
#endif
#ifdef stat_ions
     write(6,*) '     stationary ions'
#endif
#ifdef subtract_average_rho
     write(6,*) '     subtract average charge'
#endif
#ifdef electrostatic
     write(6,*) '     electrostatic simulation'
#endif
#ifdef relativistic
     write(6,*) '     relativistic simulation'
#endif
     write(6,*) '  boundary conditions:'
#if (boundary_condition == periodic)
     write(6,*) '       periodic'
#elif (boundary_condition == gem_reconnection)
     write(6,*) '       gem_reconnection'
#elif (boundary_condition == gem_3d)
     write(6,*) '       gem_3d'
#else
  we never should get here
#endif
     write(6,*) '  multigrid parameters:'
     write(6,*) '       nu1 = ',nu1
     write(6,*) '       nu2 = ',nu2
     write(6,*) '       nu3 = ',nu3
     write(6,*) '       eps = ',eps
#ifdef alt_restrict_prolong
     write(6,*) '       alt_restrict_prolong'
#endif
#ifdef coarsegrid_sor
     write(6,*) '       SOR at coarsest grid'
#endif
#ifdef miniter
     write(6,*) '       miniter =',miniter
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
     write(6,*) '       skip_poisson =', skip_poisson
     write(6,*) '  default precision       = ',drk
     write(6,*) '  precision for particles = ',prk
     write(6,*) '********** parameters **********'
  endif

! read start-up file
  call input(filecount,.true.)
  call partsort(rvi,np_i)
  call partsort(rve,np_e)
! print load balance statistics
  call load_balance()
  t_diag=t+diagout
! initialize e2,b2
  do z=1,nz; do y=1,ny; do x=1,nx
    b2x(x,y,z)=b1x(x,y,z); b2y(x,y,z)=b1y(x,y,z); b2z(x,y,z)=b1z(x,y,z)
    e2x(x,y,z)=e1x(x,y,z); e2y(x,y,z)=e1y(x,y,z); e2z(x,y,z)=e1z(x,y,z)
  enddo; enddo; enddo
! initialize rho,j
  call calc_rho(rvi,np_i,real(+1.,drk),.true.)
  call calc_rho(rve,np_e,real(-1.,drk),.false.)
#ifdef subtract_average_rho
  call subtract_average(rho)
#endif
  call smooth_rho(rho,rho)
#ifndef electrostatic
  call calc_j(rvi,np_i,real(+1.,drk),.true.)
  call calc_j(rve,np_e,real(-1.,drk),.false.)
  call smooth_j(jx,jy,jz,jx,jy,jz)
#endif
! fix divergence of b and e
  call divergence(e1x,e1y,e1z,b1x,b1y,b1z,.true.)
!  call divergence(e2x,e2y,e2z,b2x,b2y,b2z,.true.)
! NEW: B must be smoothed initially to keep j = curl(B).  If not done
! a light wave is launched in the system.
  call smooth_b(b1x,b1y,b1z,b1x,b1y,b1z)
  call smooth_b(b2x,b2y,b2z,b2x,b2y,b2z)
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
#ifndef electrostatic
    do i=1,substeps 
      call stepfield()
    enddo
#endif
#ifndef stat_ions
    call stepx(rvi,np_i)
    if(mod(n,10)==0) call partsort(rvi,np_i)
#endif
    call stepx(rve,np_e)
    if(mod(n,10)==0) call partsort(rve,np_e)
    if (mod(n,skip_poisson)==0) then
      call calc_rho(rvi,np_i,real(+1.,drk),.true.)
      call calc_rho(rve,np_e,real(-1.,drk),.false.)
#ifdef subtract_average_rho
      call subtract_average(rho)
#endif
      call smooth_rho(rho,rho)
      call divergence(e1x,e1y,e1z,b1x,b1y,b1z,.false.)
!      call divergence(e2x,e2y,e2z,b2x,b2y,b2z,.false.)
    endif
#ifndef stat_ions
    call stepv(rvi,np_i,real(1.,drk))
    call stepx(rvi,np_i)
#endif
    call stepv(rve,np_e,real(-m_e,drk))
    call stepx(rve,np_e)
#ifndef electrostatic
    call calc_j(rvi,np_i,real(+1.,drk),.true.)
    call calc_j(rve,np_e,real(-1.,drk),.false.)
    call smooth_j(jx,jy,jz,jx,jy,jz)
    do i=1,substeps 
     call stepfield()
   enddo
#endif
    t=t+(dt)
#define movie_write
#include movie_header
#undef movie_write
    if (t > t_diag) then
! print load balance statistics
      call load_balance()
      t_diag=t+diagout
    endif
! check, if time limit reached
!    presenttime=MPI_Wtime()
!    call MPI_Bcast(presenttime,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
!    if (presenttime .ge. stoptime) then
!      if (myproc == 0) write (6,*) 'time limit reached, save data'
!      terminate=.true.
!    endif
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
!  if (myproc == 0) then
#define movie_close
#include movie_header
#undef movie_close
!  endif
  call MPI_Finalize(mpi_err)
End Program p3d

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
