!******************************************************************************
!        3D particle code: double electron sheet initialization
!                           Andreas Zeiler, 2000
!
!                       LATEST CHANGE: July 17, 2001
!
!******************************************************************************

! initelectron and initelectron3: current along z
! initelectron2 and initelectron4: current along x
! initelectron and initelectron2: no ions
! initelectron3 and initelectron4: homogeneous ions

#include "param"
#define initelectron 1
#define initelectron2 2
#define initelectron3 3
#define initelectron4 4

#ifndef b0
#define b0 0.
#endif
#ifndef psi0
#define psi0 0.
#endif
#ifndef w0
#define w0 1.
#endif
#ifndef T_e
#define T_e 0.
#endif
#ifndef T_i
#define T_i 0.
#endif

Subroutine init_electron_all()
  Use pe_env
  Use partfield
  Implicit None

  Real(kind=drk) gx,gy,gz,yy,b,n,ntotal,ey_bx,rand_num
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez),&
                               pi2=6.2831853, pi4=2*pi2,ddy=dy*0.01
  Integer i,x,y,z,pe,sign,nn

  gx(x)=((x-0.5)+my_pex*nx)*dx
  gy(y)=((y-0.5)+my_pey*ny)*dy
  gz(z)=((z-0.5)+my_pez*nz)*dz

! shape of the sheet
  b(yy)=tanh((yy-ly*.25)/w0)-tanh((yy-ly*.75)/w0)-1
  n(yy)=1+(b(yy+ddy)**2-2*b(yy)**2+b(yy-ddy)**2)/(ddy**2*2*c_2)

  call init_pe_env()
  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** init: compiled for different number of PEs *****'
    call exitallpes()
  endif

! write parameters to log-file
  if (myproc==0) then
#if (init_scheme == initelectron)
    write(6,*) '******** initelectrons *********'
#elif (init_scheme == initelectron2)
    write(6,*) '******** initelectrons2 ********'
#elif (init_scheme == initelectron3)
    write(6,*) '******** initelectrons3 ********'
#elif (init_scheme == initelectron4)
    write(6,*) '******** initelectrons4 ********'
#endif
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
    write(6,*) '     c_2 = ',c_2
    write(6,*) '     part./gridp. = ',ppg
    write(6,*) '     T_e = ',T_e
    write(6,*) '     T_i = ',T_i
    write(6,*) '     b0 = ',b0
    write(6,*) '     w0 = ',w0
    write(6,*) '     psi0 = ',psi0
#ifdef relativistic
    write(6,*) '     relativistic initialization'
#endif
    write(6,*) '********** parameters **********'
  endif

! set electric and magnetic field
  do z=1,nz; do y=1,ny; do x=1,nx
#if (init_scheme == initelectron || init_scheme == initelectron3)
    b1x(x,y,z)=b(gy(y)); b1y(x,y,z)=0.; b1z(x,y,z)=b0
#elif (init_scheme == initelectron2 || init_scheme == initelectron4)
    b1x(x,y,z)=-b0; b1y(x,y,z)=0.; b1z(x,y,z)=b(gy(y))
#endif
    e1x(x,y,z)=0.; e1y(x,y,z)=0.; e1z(x,y,z)=0.
  enddo; enddo; enddo 
  np_i=0; np_e=0
! load electrons
  do i=1,nx*ny*nz*ppg*20
    np_e=np_e+1
    if (np_e .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
    endif
    call location(rve(1,np_e),rve(2,np_e),rve(3,np_e))
!        accept or reject particle according to local density
    if (rand_num()*20 .ge. n(real(rve(2,np_e),drk))) then
      np_e=np_e-1
    else
      call maxwellian(real(T_e,drk), real(m_e,drk), &
                      rve(4,np_e), rve(5,np_e), rve(6,np_e))
!        add initial electron velocity
      ey_bx=-(b(rve(2,np_e)+ddy)-b(rve(2,np_e)-ddy))/(2*ddy)
#if (init_scheme == initelectron || init_scheme == initelectron3)
    call veloadd(rve(4,np_e),rve(5,np_e),rve(6,np_e), &
                      real(0.,drk), real(0.,drk), -ey_bx)
#elif (init_scheme == initelectron2 || init_scheme == initelectron4)
    call veloadd(rve(4,np_e),rve(5,np_e),rve(6,np_e), &
                      ey_bx, real(0.,drk), real(0.,drk))
#endif
    endif
  enddo
  n_avg=ppg
  call redistribute(rve,np_e)
#if (init_scheme == initelectron3 || init_scheme == initelectron4)
! load ions
  ntotal=np_e
  ntotal=totalsum(ntotal)
  np_i=int(ntotal/(pex*pey*pez)+0.5)
  if (myproc == 0) np_i=int(ntotal-np_i*(pex*pey*pez-1)+0.5)
  if (np_i .gt. maxparticles) then
    write(6,*) '***** init: particle buffer overflow *****'
    call exitallpes()
  endif
  do i=1,np_i
    call location(rvi(1,i),rvi(2,i),rvi(3,i))
    call maxwellian(real(T_i,drk), real(1.,drk), rvi(4,i), rvi(5,i), rvi(6,i))
  enddo
  call redistribute(rvi,np_i)
#endif
! calculate electric field
#if (init_scheme == initelectron || init_scheme == initelectron2)
  call calc_rho(rve,np_e,real(-1.,drk),.true.)
#elif (init_scheme == initelectron3 || init_scheme == initelectron4)
  call calc_rho(rvi,np_i,real(+1.,drk),.true.)
  call calc_rho(rve,np_e,real(-1.,drk),.false.)
#endif
  call subtract_average(rho)
  call smooth_rho(rho,rho)
  call divergence(e1x,e1y,e1z,b1x,b1y,b1z,.false.)
#if (init_scheme == initelectron || init_scheme == initelectron3)
! add perturbation to magnetic field to start reconnection
  do z=1,nz; do y=1,ny; do x=1,nx
    b1x(x,y,z)=b1x(x,y,z)-psi0*pi4/ly*sin(pi2*gx(x)/lx)*sin(pi4*gy(y)/ly)
    b1y(x,y,z)=b1y(x,y,z)+psi0*pi2/lx*cos(pi2*gx(x)/lx)*(1-cos(pi4*gy(y)/ly))
  enddo; enddo; enddo 
#endif
! write output file
  call output(0,.true.)
End Subroutine init_electron_all
