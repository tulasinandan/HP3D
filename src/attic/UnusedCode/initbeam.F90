!******************************************************************************
!             3D particle code: double beam initialization
!                           Andreas Zeiler, 2001
!
!                       LATEST CHANGE: July 16, 2001
!
!******************************************************************************

#include "param"
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
#ifndef v0
#define v0 0.
#endif

Subroutine init_beam()
  Use pe_env
  Use partfield
  Implicit None

  Real(kind=drk) gx,gy,gz,yy,b,n,rand_num,sech
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez)
  Integer i,x,y,z,pe,sign

  gx(x)=((x-0.5)+my_pex*nx)*dx
  gy(y)=((y-0.5)+my_pey*ny)*dy
  gz(z)=((z-0.5)+my_pez*nz)*dz
  sech(yy)=1/cosh(yy)

! shape of the sheet
  n(yy)=(sech(real(((yy-ly*.25)/w0),drk))**2+ &
         sech(real(((yy-ly*.75)/w0),drk))**2)
  b(yy)=sqrt(2*(1-n(yy))*(T_i+T_e)+1)

  call init_pe_env()
  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** init: compiled for different number of PEs *****'
    call exitallpes()
  endif

! write parameters to log-file
  if (myproc==0) then
     write(6,*) '*********** initbeam ***********'
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
     write(6,*) '     m_e = ',m_e
     write(6,*) '     part./gridp. = ',ppg
     write(6,*) '     T_i = ',T_i
     write(6,*) '     T_e = ',T_e
#ifndef no_background_density
     write(6,*) '     n_0 = ',n_0
#else
     write(6,*) '     n_0 = 0'
#endif
     write(6,*) '     v0 = ',v0
     write(6,*) '     w0 = ',w0
#ifdef relativistic
     write(6,*) '     relativistic initialization'
#endif
     write(6,*) '********** parameters **********'
  endif

! set electric and magnetic field
  do z=1,nz; do y=1,ny; do x=1,nx
    b1x(x,y,z)=b(gy(y)); b1y(x,y,z)=0.; b1z(x,y,z)=0.
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
! load ions in beams
  do i=1,nx*ny*nz*ppg*2/n_0
    np_i=np_i+1
    if (np_i .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
    endif
!        accept or reject particle according to local density
    call location(rvi(1,np_i),rvi(2,np_i),rvi(3,np_i))
    if (rand_num()*2 .ge. n(real(rvi(2,np_i),drk))) then
      np_i=np_i-1
    else
      call maxwellian(real(T_i,drk), real(1.,drk), &
                      rvi(4,np_i), rvi(5,np_i), rvi(6,np_i))
      sign=1; if (rvi(2,np_i) .gt. ly/2) sign=-1
      call veloadd(rvi(4,np_i),rvi(5,np_i),rvi(6,np_i), &
                   real(v0,drk), real(0.,drk), real(0.,drk) )
    endif
  enddo
! load electrons in beams
  do i=1,nx*ny*nz*ppg*2/n_0
    np_e=np_e+1
    if (np_e .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
    endif
!        accept or reject particle according to local density
    call location(rve(1,np_e),rve(2,np_e),rve(3,np_e))
    if (rand_num()*2 .ge. n(real(rve(2,np_e),drk))) then
      np_e=np_e-1
    else
      call maxwellian(real(T_e,drk), real(m_e,drk), &
                      rve(4,np_e), rve(5,np_e), rve(6,np_e))
      sign=1; if (rve(2,np_e) .gt. ly/2) sign=-1
      call veloadd(rve(4,np_e),rve(5,np_e),rve(6,np_e), &
                     real(v0,drk), real(0.,drk), real(0.,drk) )
    endif
  enddo
  n_avg=ppg/n_0
  call redistribute(rvi,np_i)
  call redistribute(rve,np_e)
  call output(0,.true.)
End Subroutine init_beam
