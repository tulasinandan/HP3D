!******************************************************************************
!            3D particle code: single harris sheet initialization 
!                           Andreas Zeiler, 1999
!
!                        LATEST CHANGE: July 16, 2001
!
!******************************************************************************

#include "param"
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

Subroutine init_test()
  Use pe_env
  Use partfield
  Implicit None

  Real(kind=drk) :: gx,gy,gz,yy,b,n,rand_num,sech,grad_n,max_pplane,&
	denom,grad_pe,grad_pi
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez),&
                               pi2=6.2831853, pi1=pi2/2,pi4 =2*pi2
  Integer i,x,y,z,pe,nsup
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: curl_x,curl_y,curl_z,square_b,jdotb
  Real(kind=drk) :: realx,realy,realz
  Integer :: invgx,invgy,invgz

  gx(x)=((x-0.5)+my_pex*nx)*dx
  gy(y)=((y-0.5)+my_pey*ny)*dy
  gz(z)=((z-0.5)+my_pez*nz)*dz  
  sech(yy)=1/cosh(yy)

! the grid point onto which realx maps
 invgx(realx) = int(mod(realx/dx,real(nx,drk)))+1
 invgy(realy) = int(mod(realy/dy,real(ny,drk)))+1
 invgz(realz) = int(mod(realz/dz,real(nz,drk)))+1

! shape of the sheet
  b(yy)=tanh((yy-ly*.25)/w0)-tanh((yy-ly*.75)/w0)-1.
  n(yy)= -0.45*(tanh((yy-ly*.25)/w0)-tanh((yy-ly*.75)/w0))+1.
  grad_n(yy) = -0.45/w0*((sech((yy-ly*.25)/w0))**2 - &
	(sech((yy-ly*.75)/w0))**2)
  nsup =int(1./(T_i+T_e) + 1.)

  call init_pe_env()
  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** init: compiled for different number of PEs *****'
    call exitallpes()
  endif

! write parameters to log-file
  if (myproc==0) then
     write(6,*) '**********  inittest   **********'
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
!#ifndef no_background_density
     write(6,*) '     n_0 = ',n_0
!#else
!     write(6,*) '     n_0 = 0'
!#endif
     write(6,*) '     b0 = ',b0
     write(6,*) '     w0 = ',w0
     write(6,*) '     psi0 = ',psi0
#ifdef relativistic
     write(6,*) '     relativistic initialization'
#endif
     write(6,*) '********** parameters **********'
  endif

  do z=1,nz; do y=1,ny; do x=1,nx
    b1x(x,y,z)=b(gy(y)); b1y(x,y,z)=0.
  enddo; enddo; enddo

! THIS IS NOT REALLY E1X.  I AM USING THE VARIABLE SPACE TO DEFINE
! THE INPLANE PRESSURE
  do z = 1,nz; do y = 1,ny; do x = 1,nx
    e1x(x,y,z) = n(gy(y))*(T_i+T_e) + 0.5*(b1x(x,y,z)**2+b1y(x,y,z)**2)
  enddo; enddo; enddo

! calculate global minimum and maximum
  max_pplane = -1e20
  do z=1,nz; do y=1,ny; do x=1,nx
      max_pplane=max(max_pplane,e1x(x,y,z))
  enddo; enddo; enddo
  max_pplane=totalmax(max_pplane)

  do z=1,nz; do y=1,ny; do x=1,nx
    b1z(x,y,z)=sqrt(b0**2 + 2.*(max_pplane - e1x(x,y,z)))
    square_b(x,y,z) = b1x(x,y,z)**2 + b1y(x,y,z)**2 + b1z(x,y,z)**2
  enddo; enddo; enddo

! add perturbation to magnetic field to start reconnection
  do z=1,nz; do y=1,ny; do x=1,nx
    b1x(x,y,z)=b1x(x,y,z)-psi0*pi4/ly*sin(pi2*gx(x)/lx)*sin(pi4*gy(y)/ly)
    b1y(x,y,z)=b1y(x,y,z)+psi0*pi2/lx*cos(pi2*gx(x)/lx)*(1-cos(pi4*gy(y)/ly))
  enddo; enddo; enddo

  do z=1,nz; do y=1,ny; do x=1,nx
    square_b(x,y,z) = b1x(x,y,z)**2 + b1y(x,y,z)**2 + b1z(x,y,z)**2
  enddo; enddo; enddo

!this updates the ghost_cells x=0,x=nx+1,...
  call bound_b(b1x,b1y,b1z)
  do z=1,nz; do y=1,ny; do x=1,nx
    curl_x(x,y,z)=(b1z(x,y+1,z)-b1z(x,y-1,z))/(2*dy)- &
                  (b1y(x,y,z+1)-b1y(x,y,z-1))/(2*dz) 
    curl_y(x,y,z)=(b1x(x,y,z+1)-b1x(x,y,z-1))/(2*dz)- &
                  (b1z(x+1,y,z)-b1z(x-1,y,z))/(2*dx) 
    curl_z(x,y,z)=(b1y(x+1,y,z)-b1y(x-1,y,z))/(2*dx)- &
                  (b1x(x,y+1,z)-b1x(x,y-1,z))/(2*dy)
  enddo; enddo; enddo

! Note that e1x is redfined here.
  do z=1,nz; do y=1,ny; do x=1,nx
    jdotb(x,y,z)=curl_x(x,y,z)*b1x(x,y,z)+curl_y(x,y,z)*b1y(x,y,z) +&
	curl_z(x,y,z)*b1z(x,y,z)
    e1x(x,y,z)=0.
    e1y(x,y,z)=0.
    e1z(x,y,z)=0.
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
! accept or reject particle according to local density
    call location(rvi(1,np_i),rvi(2,np_i),rvi(3,np_i))
rvi(1,np_i) = nx*dx*(my_pex + rvi(1,np_i)/lx)
rvi(2,np_i) = ny*dy*(my_pey + rvi(2,np_i)/ly)
rvi(3,np_i) = nz*dz*(my_pez + rvi(3,np_i)/lz)
    if (rand_num()*nsup .ge. n(real(rvi(2,np_i),drk))) then
      np_i=np_i-1
    else
      call maxwellian(real(T_i,drk), real(1.,drk),&
	rvi(4,np_i),rvi(5,np_i),rvi(6,np_i))
      denom = n(real(rvi(2,np_i),drk))*square_b(invgx(real(rvi(1,np_i),drk)),&
	invgy(real(rvi(2,np_i),drk)),invgz(real(rvi(3,np_i),drk)))
      grad_pi = T_i*grad_n(real(rvi(2,np_i),drk))
      call veloadd(rvi(4,np_i),rvi(5,np_i),rvi(6,np_i), &
	(-grad_pi*b1z(invgx(real(rvi(1,np_i),drk)),invgy(real(rvi(2,np_i),drk)),&
	invgz(real(rvi(3,np_i),drk))))/denom, real(0.,drk),&
	(grad_pi*b1x(invgx(real(rvi(1,np_i),drk)),invgy(real(rvi(2,np_i),drk)),&
	invgz(real(rvi(3,np_i),drk))))/denom)
    endif 
  enddo
! load electrons in current sheet
  do i=1,nx*ny*nz*ppg*nsup/n_0
    np_e=np_e+1
    if (np_e .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
    endif
! accept or reject particle according to local density
    call location(rve(1,np_e),rve(2,np_e),rve(3,np_e))
rve(1,np_e) = nx*dx*(my_pex + rve(1,np_e)/lx)
rve(2,np_e) = ny*dy*(my_pey + rve(2,np_e)/ly)
rve(3,np_e) = nz*dz*(my_pez + rve(3,np_e)/lz)
    if (rand_num()*nsup .ge. n(real(rve(2,np_e),drk))) then
      np_e=np_e-1
    else
      call maxwellian(real(T_e,drk), real(m_e,drk), &
                      rve(4,np_e), rve(5,np_e), rve(6,np_e))
      denom = n(real(rve(2,np_e),drk))*square_b(invgx(real(rve(1,np_e),drk)),&
	invgy(real(rve(2,np_e),drk)),invgz(real(rve(3,np_e),drk)))
      grad_pe = T_e*grad_n(real(rve(2,np_e),drk))
      call veloadd(rve(4,np_e),rve(5,np_e),rve(6,np_e), &
	(grad_pe*b1z(invgx(real(rve(1,np_e),drk)),invgy(real(rve(2,np_e),drk)),&
	invgz(real(rve(3,np_e),drk)))-b1x(invgx(real(rve(1,np_e),drk)),&
	invgy(real(rve(2,np_e),drk)),invgz(real(rve(3,np_e),drk)))*&
	jdotb(invgx(real(rve(1,np_e),drk)),invgy(real(rve(2,np_e),drk)),&
	invgz(real(rve(3,np_e),drk))))/denom, &
	-b1y(invgx(real(rve(1,np_e),drk)),&
	invgy(real(rve(2,np_e),drk)),invgz(real(rve(3,np_e),drk)))*&
	jdotb(invgx(real(rve(1,np_e),drk)),invgy(real(rve(2,np_e),drk)),&
	invgz(real(rve(3,np_e),drk)))/denom,&
	(-grad_pe*b1x(invgx(real(rve(1,np_e),drk)),invgy(real(rve(2,np_e),drk)),&
	invgz(real(rve(3,np_e),drk)))-b1z(invgx(real(rve(1,np_e),drk)),&
	invgy(real(rve(2,np_e),drk)),invgz(real(rve(3,np_e),drk)))*&
	jdotb(invgx(real(rve(1,np_e),drk)),invgy(real(rve(2,np_e),drk)),&
	invgz(real(rve(3,np_e),drk))))/denom)
    endif
  enddo
  n_avg=ppg/n_0
  call redistribute(rvi,np_i)
  call redistribute(rve,np_e)

  call output(0,.true.)
End Subroutine init_test

