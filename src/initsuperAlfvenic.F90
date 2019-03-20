!================================================================================
!
!
!
!  Initialize a superalfvenic flow, similar to SW, the hard coded value is 
!  8.0*Va in the x direction. The magnetic field direction is decided by the
!  values guven in the param file
!
!		Tulasi Nandan Parashar 09/11/2013
!
!
!
!================================================================================

#include "param"
#ifndef b0x
#define b0x 0.
#endif
#ifndef b0y
#define b0y 0.
#endif
#ifndef b0z
#define b0z 0.
#endif
#ifndef psi0
#define psi0 0.
#endif
#ifndef w0
#define w0 1.
#endif
#ifndef T_i
#define T_i 0.
#endif
#ifndef T_e
#define T_e 0.
#endif
#ifndef n_0
#define n_0 1.
#define no_background_density
#endif

Subroutine init_superAlfvenic()
  Use pe_env
  Use partfield
  Implicit None

  Real(kind=drk) gx,gy,gz,xx,yy,zz,b,n,rand_num,sech,ppert,npert,va
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez),&
                               pi2=6.2831853, pi4=2*pi2
  Integer i,x,y,z,pe,sign,nsup
! Tulas's parameters for wave 09/20/2006
  Real(kind=drk) :: sum_ncell, ntot_p
  Integer :: mpi_err
  Integer :: part_counter
  gx(x)=((x-0.5)+my_pex*nx)*dx
  gy(y)=((y-0.5)+my_pey*ny)*dy
  gz(z)=((z-0.5)+my_pez*nz)*dz
  sech(yy)=1/cosh(yy)

! ppert(xx,yy)=0.
  ppert(xx,yy)=0.

  npert(xx,yy)= ppert(xx,yy)/T_i
  nsup=int(1./(T_i+T_e)+5.)

  call init_pe_env()
  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** init-hybrid: compiled for different number of PEs *****'
    call exitallpes()
  endif

!write parameters to log-file
  if (myproc==0) then
    write(6,*) '********** initnone **********'
    write(6,*) '******* (hybrid version) *******'
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
    write(6,*) '     part./gridp. = ',ppg
    write(6,*) '     T_i = ',T_i
    write(6,*) '     T_e = ',T_e
#ifndef no_background_density
    write(6,*) '     n_0 = ',n_0
#else
    write(6,*) '     n_0 = 0'
#endif
    write(6,*) '     b0 = ',b0x, b0y, b0z
    write(6,*) '     w0 = ',w0
    write(6,*) '     psi0 = ',psi0
    write(6,*) '********** parameters **********'
  endif

! set magnetic field
  do z=1,nz; do y=1,ny; do x=1,nx
    b1x(x,y,z)=b0x; b1y(x,y,z)=b0y; b1z(x,y,z)=b0z
  enddo; enddo; enddo 

! Find the value of Alfven speed va
  va = sqrt(b0x**2+b0y**2+b0z**2)/sqrt(n_0)
!-------------------------------------------
! finding the total number of particles
! Ntot_proc=n_avg*Sum(n_cell)/(pex*pey)
!	 		Tulasi
  ntot_p=0
  do x=1,nx; do y=1,ny
    ntot_p = ntot_p + npert(gx(x),gy(y)) 
  enddo; enddo
  ntot_p = nz*ntot_p
  
  call MPI_Allreduce(ntot_p,sum_ncell,1,mpi_drk,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  n_avg=ppg/n_0 
!-------------------------------------------
  np_i=0
  part_counter=0
! load ions in current sheets
!
 inf_loop: do
 !old  do i=1,nx*ny*nz*ppg*nsup/n_0
    np_i=np_i+1
    if (np_i .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
 !old    endif
 ! This is to make sure that we get the number of particles
 ! we want. We did not believe the previous version.
    elseif (np_i .ge. ppg*nx*ny*nz+int(ntot_p*n_avg)) then
    exit inf_loop
    else
 !        accept or reject particle according to local density
    call location(rvi(1,np_i),rvi(2,np_i),rvi(3,np_i))
    if (rand_num()*nsup .ge. n_0+npert(real(rvi(1,np_i),drk),real(rvi(2,np_i),drk))) then
! checking the number of particles if loop
      np_i=np_i-1
    else
      call maxwellian(real(T_i,drk), real(1.,drk), &
                      rvi(4,np_i), rvi(5,np_i), rvi(6,np_i))

      if ((myproc .eq. 23) .and. (gx(32) .lt. real(rvi(1,np_i),drk)) .and. &
          (real(rvi(1,np_i),drk) .lt. gx(34)) .and. (gy(32) .lt. real(rvi(2,np_i),drk)) &
          .and. (real(rvi(2,np_i),drk) .lt. gy(34))) then
       part_counter=part_counter+1
      endif
      sign=1; if (rvi(2,np_i) .gt. ly/2) sign=-1
      call veloadd(rvi(4,np_i),rvi(5,np_i),rvi(6,np_i), &
!        real(-sin(pi2*rvi(2,np_i)/ly),drk), real(sin(pi2*rvi(1,np_i)/lx),drk), real(0.,drk) )
         real(8.*va,drk), real(0.,drk), real(0.,drk) )
!we add the velocities 
    endif ! end the check number of particles loop
    endif ! for the maxparticles if loop
 !old  enddo
    enddo inf_loop
!-------------------------------------------

   n_avg=ppg/n_0  
  call redistribute(rvi,np_i)

! calculate electron pressure
    call calc_rho(rvi,np_i,real(+1.,drk),.true.)
  do z=1,nz; do y=1,ny; do x=1,nx
#ifdef barotropic
    pe1(x,y,z)=T_e*rho(x,y,z)**gam                
#else
    pe1(x,y,z)=rho(x,y,z)*T_e
#endif
  enddo; enddo; enddo 
! write start-up file
  call output(0,.true.)
End Subroutine init_superAlfvenic
