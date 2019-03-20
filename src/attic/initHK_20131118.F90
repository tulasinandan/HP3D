!================================================================================
!
!
!
!   This is the Alfven wave initialization subroutine for the hybrid code.
!
!
!		Tulasi Nandan Parashar 09/20/2006
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
#ifndef b0y
#define b0y 0.
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

Subroutine init_HK()
  Use pe_env
  Use partfield
  Implicit None

  Real(kind=drk) gx,gy,gz,yy,b,n,rand_num,sech
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez),&
                               pi2=6.2831853, pi4=2*pi2
  Integer i,x,y,z,pe,sign,nsup
! Tulasi's parameters for Alfven wave 09/20/2006
  Real(kind=drk) cak, speedx, speedy, speedz
  Real(kind=drk), Parameter :: k_wave=nummodes*pi2/lx

  gx(x)=((x-0.5)+my_pex*nx)*dx
  gy(y)=((y-0.5)+my_pey*ny)*dy
  gz(z)=((z-0.5)+my_pez*nz)*dz
  sech(yy)=1/cosh(yy)

! Initial condition
  !Propagation theta
  b(yy)=waveamp*cos(k_wave*yy)
  n(yy)=0.
  cak(yy) = cos(atan(b0y/b0x))*sqrt((b0x**2+b0y**2+b0z**2)/(n(yy)+n_0))
! speedx(yy)=waveamp*cos(k_wave*yy)*cak(yy)*b0z/((1-(T_i+T_e)/cak(yy)**2)*b0x)
  speedx(yy)=0.
! speedy(yy)=waveamp*sin(k_wave*yy)
! speedz(yy)=-waveamp*cos(k_wave*yy)*cak(yy)
  speedy(yy)=-waveamp*cos(k_wave*yy)*cak(yy)
  speedz(yy)=waveamp*sin(k_wave*yy)
! speed(yy)=5
  nsup=int(1./(T_i+T_e)+1.)

  call init_pe_env()
  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** init-hybrid: compiled for different number of PEs *****'
    call exitallpes()
  endif

!'write parameters to log-file
  if (myproc==0) then
    write(6,*) '********** initalfven **********'
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
!   b1x(x,y,z)=b0x; b1y(x,y,z)=b0y; b1z(x,y,z)=b(gx(x))+b0z
    b1x(x,y,z)=b0x; b1y(x,y,z)=b(gx(x))+b0y; b1z(x,y,z)=b0z
  enddo; enddo; enddo 
  np_i=0
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
    elseif (np_i .gt. (nx*ny*nz*ppg)) then
    exit inf_loop
    else
 !        accept or reject particle according to local density
    call location(rvi(1,np_i),rvi(2,np_i),rvi(3,np_i))
      call maxwellian(real(T_i,drk), real(1.,drk), &
                      rvi(4,np_i), rvi(5,np_i), rvi(6,np_i))
!   rvi(4,np_i)=0.; rvi(5,np_i)=0.; rvi(6,np_i)=0.
      sign=1; if (rvi(2,np_i) .gt. ly/2) sign=-1
      call veloadd(rvi(4,np_i),rvi(5,np_i),rvi(6,np_i), &
                   real(speedx(rvi(1,np_i)),drk), &
                   real(speedy(rvi(2,np_i)),drk), &
                   real(speedz(rvi(1,np_i)),drk) )

    endif ! for the maxparticles if loop
 !old  enddo
    enddo inf_loop

   n_avg=ppg/n_0  
  call redistribute(rvi,np_i)
! add perturbation to magnetic field to start reconnection
! do z=1,nz; do y=1,ny; do x=1,nx
!   b1x(x,y,z)=b1x(x,y,z)-psi0*pi4/ly*sin(pi2*gx(x)/lx)*sin(pi4*gy(y)/ly)
!   b1y(x,y,z)=b1y(x,y,z)+psi0*pi2/lx*cos(pi2*gx(x)/lx)*(1-cos(pi4*gy(y)/ly))
! enddo; enddo; enddo 
! calculate electron pressure
  do z=1,nz; do y=1,ny; do x=1,nx
#ifdef no_background_density
    pe1(x,y,z)=n(gy(y))*T_e
#else
    pe1(x,y,z)=(n(gy(y))+n_0)*T_e
#endif
  enddo; enddo; enddo 
! write start-up file
  call output(0,.true.)
End Subroutine init_HK
