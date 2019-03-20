!******************************************************************************
!              3D particle code: main module, hybrid version
!                           Andreas Zeiler, 2001
!
!                       LATEST CHANGE: July 20, 2001
!
!******************************************************************************

#include "param"
#ifndef hybrid
  switch hybrid must be set for hybrid version
#endif

#ifndef D_b
#define D_b 0.
#endif
#ifndef D_pe
#define D_pe 0.
#endif
#ifndef GS_LEX
#ifndef GS_RB
#ifndef JACOBI
#define JACOBI
#endif
#endif
#endif
#ifndef take_dump
#define take_dump 28800
#endif

! assign unique number to each set of boundary conditions
#define periodic 0

Program p3d_hybrid
  Use pe_env
  Use partfield
  Use timer
#ifdef EXPANSION
  Use expansion_params
#endif
  Implicit none

  Integer :: n,nsteps,i,filecount,x,y,z
  Character*16 :: buf
#ifdef EXPANSION
  Character*16 :: junk
#endif
  Real(kind=drk) :: t_diag
  Real(kind=drk) :: stoptime, presenttime
  Real(kind=drk) :: movie_dmp_time
  Logical :: terminate=.false.
  Integer :: mpi_err
#ifdef sn6311
  Logical, External :: mpp_terminate
#endif
#define movie_variables
#include movie_header
#include movie_header_full
#undef movie_variables

  call init_pe_env()
! if (nprocs .ne. n_pes .and. myproc == 0) then
!   write(6,*) '***** p3d-hybrid: compiled for different number of PEs *****'
!   call exitallpes()
! endif

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
     write(6,*) '********** p3d-hybrid **********'
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
     write(6,*) '     gamma = ',gam
#ifdef n_base
     write(6,*) '     n_base = ',n_base
#endif
#ifdef isothermal
     write(6,*) '     isothermal simulation'
     write(6,*) '        T_e = ',T_e
#endif
     write(6,*) '     dt = ',dt
     write(6,*) '     substeps = ',substeps
#ifdef double_byte
     write(6,*) '     two-byte movie output'
#endif
     write(6,*) '  boundary conditions:'
#if (boundary_condition == periodic)
     write(6,*) '       periodic'
#else
  we never should get here
#endif
     write(6,*) '  hyperviscosities:'
     write(6,*) '       D_b  = ',D_b
     write(6,*) '       D_pe = ',D_pe
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
     write(6,*) '  default precision       = ',drk
     write(6,*) '  precision for particles = ',prk
     write(6,*) '********** parameters **********'
  endif

! set the driving fields
#ifdef DRIVE_B
  if (myproc==0) print *, 'Loading the driving magnetic field'
  call set_drive_fields()
#endif
#ifdef DRIVE_V
#endif
#ifdef ACCELERATED_EXPANSION
      open(73,file='winddata.dat',form='formatted',status='old',&
       access='sequential')
      read(73,*) junk
      if (t==0) then
      read(73,*) ti,ai,daoai,wi,dwowi,usi
      read(73,*) tf,af,daoaf,wf,dwowf,usf
      dtwd=tf-ti
      endif
#endif
!BEGIN OPEN ENERGY FILE AND WRITE HEADER
    if (myproc == 0) then
    open(unit=9,file='Energies.dat',status='unknown',form='formatted')
    write(9,'(5A85)') '# 1)t       2)TOT         3)MAG      4)MIX      5)MIY      6)MIZ      7)ELC   ',&
                      '   8)ION     9)IKX        10)IKY     11)IKZ     12)IFL      13)IFX   ',&
                      '  14)IFY     15)IFZ       16)EPR     17)IPR     18)IPX      19)IPY   ',&
                      '  20)IPZ     21)VIP       22)VIX     23)VIY     24)VIZ      25)BIP   ',&
                      '  26)BIX     27)BIY       28)BIZ     29)BNU     30)ETA'
    endif
! END WRITE HEADER
! read start-up file
  call input(filecount,.true.)
  call partsort(rvi,np_i)
! call partsort(rv1,np_1)
! print load balance statistics
  call load_balance()
  t_diag=t+diagout
#ifdef EXPANSION
  call expansion()
#endif
  call divergence(b1x,b1y,b1z)
! initialize b2
! and pe2 - TNP 09/09/09
  do z=1,nz; do y=1,ny; do x=1,nx
    b2x(x,y,z)=b1x(x,y,z); b2y(x,y,z)=b1y(x,y,z); b2z(x,y,z)=b1z(x,y,z)
    pe2(x,y,z)=pe1(x,y,z)
  enddo; enddo; enddo
! initialize rho (particle density), j (ion current density)
  call calc_rho(rvi,np_i,real(+1.,drk),.true.)
  call smooth_rho(rho,rho)
  call add_base_density()
  call calc_j(rvi,np_i,real(+1.,drk),.true.)
  call smooth_j(jx,jy,jz,jx,jy,jz)
  call calc_je(b1x,b1y,b1z,jx,jy,jz,jex,jey,jez,jtotx,jtoty,jtotz)
  call calc_p(rvi,np_i,real(+1.,drk),real(+1.,drk),rho,jx,jy,jz,.true.)
! intialize movie files, write first frames
#define movie_open
#include movie_header
#include movie_header_full
#undef movie_open
  if (t .eq. 0.) then
#define movie_write
#include movie_header
#include movie_header_full
#undef movie_write
  endif
! initialize dump movie timer
  movie_dmp_time=0.
  do while (abs(t) > abs(movie_dmp_time))
    movie_dmp_time = movie_dmp_time+sign(take_dump,dt)
  enddo
  if (myproc == 0) write(6,*) 'Initial movie_dmp_time = ',movie_dmp_time


! main loop
  call init_timer()
  do n=1,nsteps
#ifdef EXPANSION
    call expansion()
    if (myproc==0) write(6,*) 'Expansion ww, dwow, aa, daoa, usw', ww,dwow,aa,daoa,usw
#endif
    call energy()
    do i=1,substeps 
!
!call stephybrid all the times for debugging.
!It actually solved the problem :)
!		Tulasi  10/19/2006
!
!     call stephybrid(i .lt. substeps)
      call stephybrid(.false.)
    enddo
    call stepx(rvi,np_i)
!   call stepx(rv1,np_1)
!   if (myproc==0) write(6,*) '#$%^#$%^#$%^#$%^#$%^ DONE INITIAL CALCULATIONS #$%^#$%^#$%^#$%^#$%^#$%^'
    if(mod(n,5)==0) call partsort(rvi,np_i)
!   if(mod(n,5)==0) call partsort(rv1,np_1)
    call calc_rho(rvi,np_i,real(+1.,drk),.true.)
    call smooth_rho(rho,rho)
    call add_base_density()
    call calc_j(rvi,np_i,real(+1.,drk),.true.)
    call smooth_j(jx,jy,jz,jx,jy,jz)
    call calc_p(rvi,np_i,real(+1.,drk),real(+1.,drk),rho,jx,jy,jz,.true.)
    call calc_e()
#ifdef EXPANSION
    call expansion()
#endif
    call stepv(rvi,np_i,real(1.,drk))
!   call stepv(rv1,np_1,real(1.,drk))
    call stepx(rvi,np_i)
!   call stepx(rv1,np_1)
    call calc_rho(rvi,np_i,real(+1.,drk),.true.)
    call smooth_rho(rho,rho)
    call add_base_density()
    call calc_j(rvi,np_i,real(+1.,drk),.true.)
    call smooth_j(jx,jy,jz,jx,jy,jz)
    call calc_p(rvi,np_i,real(+1.,drk),real(+1.,drk),rho,jx,jy,jz,.true.)
    do i=1,substeps 
!     call stephybrid(i .lt. substeps)
      call stephybrid(.false.)
    enddo
! Added a subroutine to calculate the electron current
    call calc_je(b1x,b1y,b1z,jx,jy,jz,jex,jey,jez,jtotx,jtoty,jtotz)
    t=t+dt
#define movie_write
#include movie_header
#include movie_header_full
#undef movie_write
    if (t > t_diag) then
! print load balance statistics
      call load_balance()
      t_diag=t+diagout
    endif
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
! write dump files in between to make sure that we do not lose
! a lot of computer hours if the system crashed in between a run.
  if (t .ge. movie_dmp_time) then
  filecount=filecount+1
  if (myproc == 0) write(6,*) '*****************************************'
  if (myproc == 0) write(6,*) 'Outputting dump data at time ',t                 
  if (myproc == 0) write(6,*) ' to file ',filecount
  if (myproc == 0) write(6,*) '*****************************************'
  call output(filecount,.true.)
  movie_dmp_time = movie_dmp_time + take_dump
  endif
!
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
#include movie_header_full
#undef movie_close
  endif
  call MPI_Finalize(mpi_err)
End Program p3d_hybrid

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
!#### This option for franklin
#ifdef Linux
  call getarg(num,buf)
#endif
!#### End franklin option

End Subroutine our_getarg
