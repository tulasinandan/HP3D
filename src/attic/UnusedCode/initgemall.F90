!******************************************************************************
!            3D particle code: single harris sheet initialization 
!                           Andreas Zeiler, 1999
!
!                        LATEST CHANGE: July 16, 2001
!
!******************************************************************************

#include "param"
#define initgem 1
#define initgem2 2
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

Subroutine init_gem_all()
  Use pe_env
  Use partfield
  Implicit None

  Real(kind=drk) gx,gy,gz,yy,b,n,rand_num,sech
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez),&
                               pi2=6.2831853, pi1=pi2/2
  Integer i,x,y,z,pe,nsup

  gx(x)=((x-0.5)+my_pex*nx)*dx
  gy(y)=((y-0.5)+my_pey*ny)*dy
  gz(z)=((z-0.5)+my_pez*nz)*dz
  sech(yy)=1/cosh(yy)

! shape of the sheet
  b(yy)=tanh((yy-ly*.5)/w0)
  n(yy)=(sech(real(((yy-ly*.5)/w0),drk))**2)/(2*(T_i+T_e))
  nsup=int(1./(2*(T_i+T_e))+1.)

  call init_pe_env()
  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** init: compiled for different number of PEs *****'
    call exitallpes()
  endif

! write parameters to log-file
  if (myproc==0) then
     write(6,*) '**********  initgem   **********'
     write(6,*) '********** parameters **********'
     write(6,*) '     grid points nx = ',nx
     write(6,*) '     grid points ny = ',ny
     write(6,*) '     grid points nz = ',nz
     write(6,*) '     processors pex = ',pex
     write(6,*) '     processors pey = ',pey
     write(6,*) '     edge length lx = ',lx
     write(6,*) '     edge length ly = ',ly
     write(6,*) '     edge length lz = ',lz
     write(6,*) '     m_e = ',m_e
     write(6,*) '     part./gridp. = ',ppg
     write(6,*) '     T_i = ',T_i
     write(6,*) '     T_e = ',T_e
#ifndef no_background_density
     write(6,*) '     n_0 = ',n_0
#else
     write(6,*) '     n_0 = 0'
#endif
     write(6,*) '     w0 = ',w0
     write(6,*) '     psi0 = ',psi0
#ifdef relativistic
     write(6,*) '     relativistic initialization'
#endif
     write(6,*) '********** parameters **********'
  endif

! set electric and magnetic field
  do z=1,nz; do y=1,ny; do x=1,nx
#if (init_scheme == initgem)
    b1x(x,y,z)=b(gy(y)); b1y(x,y,z)=0.; b1z(x,y,z)=0.
#elif (init_scheme == initgem2)
    b1x(x,y,z)=0.; b1y(x,y,z)=0.; b1z(x,y,z)=b(gy(y))  
#endif
    e1x(x,y,z)=0.; e1y(x,y,z)=0.; e1z(x,y,z)=0.
  enddo; enddo; enddo 
  np_i=0; np_e=0
! load particles, background density
#ifndef no_background_density
  do i=1,nx*ny*nz*ppg
    np_i=np_i+1
    if (np_i .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
    endif
    call location(rvi(1,np_i),rvi(2,np_i),rvi(3,np_i))
    call maxwellian(real(T_i,drk), real(1.,drk), &
                    rvi(4,np_i), rvi(5,np_i), rvi(6,np_i))
  enddo
  do i=1,nx*ny*nz*ppg
    np_e=np_e+1
    if (np_e .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
    endif
    call location(rve(1,np_e),rve(2,np_e),rve(3,np_e))
    call maxwellian(real(T_e,drk), real(m_e,drk), &
                    rve(4,np_e), rve(5,np_e), rve(6,np_e))
  enddo
#endif
! load ions in current sheet
  do i=1,nx*ny*nz*ppg*nsup/n_0
    np_i=np_i+1
    if (np_i .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
    endif
!        accept or reject particle according to local density
    call location(rvi(1,np_i),rvi(2,np_i),rvi(3,np_i))
    if (rand_num()*nsup .ge. n(real(rvi(2,np_i),drk))) then
      np_i=np_i-1
    else
      call maxwellian(real(T_i,drk), real(1.,drk), &
                      rvi(4,np_i), rvi(5,np_i), rvi(6,np_i))
#if (init_scheme == initgem)
      call veloadd(rvi(4,np_i),rvi(5,np_i),rvi(6,np_i), &
                   real(0.,drk), real(0.,drk), real(-2*T_i/w0,drk))
#elif (init_scheme == initgem2)
      call veloadd(rvi(4,np_i),rvi(5,np_i),rvi(6,np_i), &
                   real(2*T_i/w0,drk),real(0.,drk),real(0.,drk))
#endif
    endif 
  enddo
! load electrons in current sheet
  do i=1,nx*ny*nz*ppg*nsup/n_0
    np_e=np_e+1
    if (np_e .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
    endif
!        accept or reject particle according to local density
    call location(rve(1,np_e),rve(2,np_e),rve(3,np_e))
    if (rand_num()*nsup .ge. n(real(rve(2,np_e),drk))) then
      np_e=np_e-1
    else
      call maxwellian(real(T_e,drk), real(m_e,drk), &
                      rve(4,np_e), rve(5,np_e), rve(6,np_e))
#if (init_scheme == initgem)
      call veloadd(rve(4,np_e),rve(5,np_e),rve(6,np_e), &
                   real(0.,drk), real(0.,drk), real(2*T_e/w0,drk))
#elif (init_scheme == initgem2)
      call veloadd(rve(4,np_e),rve(5,np_e),rve(6,np_e), &
                    real(-2*T_e/w0,drk),real(0.,drk),real(0.,drk))
#endif
    endif
  enddo
  n_avg=ppg/n_0
  call redistribute(rvi,np_i)
  call redistribute(rve,np_e)
! add perturbation to magnetic field to start reconnection
  do z=1,nz; do y=1,ny; do x=1,nx
    b1x(x,y,z)=b1x(x,y,z)-psi0*pi1/ly*cos(pi2*gx(x)/lx)*cos(pi1*gy(y)/ly)
    b1y(x,y,z)=b1y(x,y,z)-psi0*pi2/lx*sin(pi2*gx(x)/lx)*sin(pi1*gy(y)/ly)
  enddo; enddo; enddo 
  call output(0,.true.)
End Subroutine init_gem_all
