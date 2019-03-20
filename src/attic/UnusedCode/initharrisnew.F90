!*****************************************************************************
!             3D particle code: double harris sheet initialization 
!                           Andreas Zeiler, 1999
!
!                      LATEST CHANGE: July 16, 2001
!
!******************************************************************************

! initharris: current along z
! initharris2: current along x

#include "param"
#define initharris 1
#define initharris2 2
#ifndef b0
#define b0 0.
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

Subroutine init_harris_all()
  Use pe_env
  Use partfield
  Implicit None

  Real(kind=drk) gx,gy,gz,yy,b,n,rand_num,sech,nsup,vid,tn
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez),&
                          pi2=6.2831853, pi4=2.*pi2
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: curl_x,curl_y,curl_z,divjxb,gridn,jdotb
  Integer i,x,y,z,pe,sign
  Real(kind=drk) :: realx,realy,realz,min_n
  Integer :: invgx,invgy,invgz,xval,yval,zval
  Real(kind=drk) :: xp,yp,zp,xm,ym,zm,inter_n,inter_jx,inter_jy,inter_jz

  gx(x)=((x-0.5)+my_pex*nx)*dx
  gy(y)=((y-0.5)+my_pey*ny)*dy
  gz(z)=((z-0.5)+my_pez*nz)*dz
  sech(yy)=1/cosh(yy)

! the _nearest_ grid point onto which realx maps
!  invgx(realx) = int(realx/dx-my_pex*nx)+1
!  invgy(realy) = int(realy/dy-my_pey*ny)+1
!  invgz(realz) = int(realz/dz-my_pez*nz)+1

! shape of the sheet
  b(yy)=tanh((yy-ly*.25)/w0)-tanh((yy-ly*.75)/w0)+&
	tanh((yy-ly*1.25)/w0)-tanh((yy+ly*.25)/w0)+1
!  n(yy)=(tanh((yy-ly*.25)/w0)-tanh((yy-ly*.75)/w0))*&
!	(2-tanh((yy-ly*.25)/w0)+tanh((yy-ly*.75)/w0))/(2*(T_i+T_e))
  n(yy) = (2-(tanh((yy-ly*.25)/w0)-tanh((yy-ly*.75)/w0)+&
	tanh((yy-ly*1.25)/w0)-tanh((yy+ly*.25)/w0)+1)**2)/(2*(T_i+T_e))
!  vid(yy) = 2.*T_i/w0*(tanh((yy-ly*.25)/w0)+tanh((yy-ly*.75)/w0))/&
!	(2-tanh((yy-ly*.25)/w0)+tanh((yy-ly*.75)/w0))
!  vid(yy)=-2*T_i/w0*(sech((yy-ly*.25)/w0)**2-sech((yy-ly*.75)/w0)**2+&
!	sech((yy-ly*1.25)/w0)**2-sech((yy+ly*.25)/w0)**2)/&
!	(1-(tanh((yy-ly*.25)/w0)-tanh((yy-ly*.75)/w0)+&
!	tanh((yy-ly*1.25)/w0)-tanh((yy+ly*.25)/w0)+1)**2)

! Anything processor dependent must be done after this next line.
  call init_pe_env()
  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** init: compiled for different number of PEs *****'
    call exitallpes()
  endif

! write parameters to log-file
  if (myproc==0) then
#if (init_scheme == initharris)
    write(6,*) '********** initharris **********'
#elif (init_scheme == initharris2)
    write(6,*) '********** initharris2 *********'
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
    write(6,*) '     part./gridp. = ',ppg
    write(6,*) '     T_i = ',T_i
    write(6,*) '     T_e = ',T_e
#ifndef no_background_density
    write(6,*) '     n_0 = ',n_0
#else
    write(6,*) '     n_0 = 0'
#endif
    write(6,*) '     b0 = ',b0
    write(6,*) '     w0 = ',w0
#if (init_scheme == initharris)
    write(6,*) '     psi0 = ',psi0
#endif
#ifdef relativistic
    write(6,*) '     relativistic initialization'
#endif
    write(6,*) '********** parameters **********'
  endif

! set electric and magnetic field (no kinetic equilibrium with guide field bz)
  do z=1,nz; do y=1,ny; do x=1,nx
#if (init_scheme == initharris)
    b1x(x,y,z)=b(gy(y)); b1y(x,y,z)=0.; b1z(x,y,z)=b0
#elif (init_scheme == initharris2)
    b1x(x,y,z)=-b0; b1y(x,y,z)=0.; b1z(x,y,z)=b(gy(y))
#endif
    e1x(x,y,z)=0.; e1y(x,y,z)=0.; e1z(x,y,z)=0.
  enddo; enddo; enddo 

tn = 2.*w0

! add perturbation to magnetic field to start reconnection
  do z=1,nz; do y=1,ny; do x=1,nx
    b1x(x,y,z)=b1x(x,y,z)-psi0*pi4/ly*sin(pi2*gx(x)/lx)*sin(pi4*gy(y)/ly)
    b1y(x,y,z)=b1y(x,y,z)+psi0*pi2/lx*cos(pi2*gx(x)/lx)*(1-cos(pi4*gy(y)/ly))
  enddo; enddo; enddo 

!this updates the ghost_cells x=0,x=nx+1,...
  call bound_b(b1x,b1y,b1z)
  do z=1,nz; do y=1,ny; do x=1,nx
    curl_x(x,y,z)=(b1z(x,y+1,z)-b1z(x,y-1,z))/(2.*dy)- &
                  (b1y(x,y,z+1)-b1y(x,y,z-1))/(2.*dz) 
    curl_y(x,y,z)=(b1x(x,y,z+1)-b1x(x,y,z-1))/(2.*dz)- &
                  (b1z(x+1,y,z)-b1z(x-1,y,z))/(2.*dx) 
    curl_z(x,y,z)=(b1y(x+1,y,z)-b1y(x-1,y,z))/(2.*dx)- &
                  (b1x(x,y+1,z)-b1x(x,y-1,z))/(2.*dy)
  enddo; enddo; enddo
  call bound_j(curl_x,curl_y,curl_z)

! Calculate div of J cross B
  do z=1,nz; do y=1,ny; do x=1,nx
    divjxb(x,y,z) = (curl_y(x+1,y,z)*b1z(x+1,y,z)-curl_y(x-1,y,z)*b1z(x-1,y,z)-&
	curl_z(x+1,y,z)*b1y(x+1,y,z)+curl_z(x-1,y,z)*b1y(x-1,y,z))/(2.*dx) + &
	(curl_z(x,y+1,z)*b1x(x,y+1,z)-curl_z(x,y-1,z)*b1x(x,y-1,z)-&
	curl_x(x,y+1,z)*b1z(x,y+1,z)+curl_x(x,y-1,z)*b1z(x,y-1,z))/(2.*dy) + &
	(curl_x(x,y,z+1)*b1y(x,y,z+1)-curl_x(x,y,z-1)*b1y(x,y,z-1)-&
	curl_y(x,y,z+1)*b1x(x,y,z+1)+curl_y(x,y,z-1)*b1x(x,y,z-1))/(2.*dz)
  enddo; enddo; enddo
  call bound_phi(0,0,divjxb)

! Put density on a grid
  do z=1,nz; do y=1,ny; do x=1,nx
    gridn(x,y,z) = n(gy(y))
  enddo; enddo; enddo
  call bound_phi(0,0,gridn)
  call poisson(gridn,1/(T_e+T_i)*divjxb,real(eps,drk))
!  call dendiv(gridn,1/(T_e+T_i)*divjxb)
  call bound_phi(0,0,gridn)

! calculate global minimum and maximum
  min_n = 1e20
  do z=1,nz; do y=1,ny; do x=1,nx
      min_n=min(min_n,gridn(x,y,z))
  enddo; enddo; enddo
  min_n=totalmin(min_n)

! Make density postive definite
  if (min_n .lt. 0) then
    min_n = abs(min_n)
    do z=1,nz; do y=1,ny; do x=1,nx
      gridn(x,y,z) = gridn(x,y,z)+min_n
    enddo; enddo; enddo
  endif
  call bound_phi(0,0,gridn)

! The integer nsup must be greater than or equal to the largest value
! of the density (excluding background density, i.e., n_0) on a
! processor.
  nsup = -1e20
  do z=1,nz; do y=1,ny; do x = 1,nx
      nsup=max(nsup,gridn(x,y,z))
  enddo; enddo; enddo

! To calculate vparallel, need jdotb
  do z=1,nz; do y=1,ny; do x=1,nx
    jdotb(x,y,z)=curl_x(x,y,z)*b1x(x,y,z)+curl_y(x,y,z)*b1y(x,y,z) +&
        curl_z(x,y,z)*b1z(x,y,z)
  enddo; enddo; enddo

  
  np_i=0; np_e=0
! load particles.  As currently written puts particles anywhere
! (i.e., not necessarily on the given processor. redistribute below
! works out the details.

! load ions in current sheets
  do i=1,nx*ny*nz*ppg*nsup/n_0
    np_i=np_i+1
    if (np_i .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
    endif
! accept or reject particle according to local density
    call location(rvi(1,np_i),rvi(2,np_i),rvi(3,np_i))

! This restricts the particle to its own processor.  More elegant
! would be to change the structure of location, but that's a little
! complicated right now.
    rvi(1,np_i) = nx*dx*(my_pex + rvi(1,np_i)/lx)
    rvi(2,np_i) = ny*dy*(my_pey + rvi(2,np_i)/ly)
    rvi(3,np_i) = nz*dz*(my_pez + rvi(3,np_i)/lz)

! Interpolate density.  First find gridpoint _below_ point.
    xval = int(rvi(1,np_i)/dx-my_pex*nx+0.5)
    yval = int(rvi(2,np_i)/dy-my_pey*ny+0.5)
    zval = int(rvi(3,np_i)/dz-my_pez*nz+0.5)
    xm=rvi(1,np_i)/dx-my_pex*nx+0.5-xval; xp=1.-xm
    ym=rvi(2,np_i)/dy-my_pey*ny+0.5-yval; yp=1.-ym
#if (nz == 1 && pez == 1)
    zm=1.; zp=0.
#else
    zm=rvi(3,np_i)/dz-my_pez*nz+0.5-zval; zp=1.-zm
#endif
    inter_n = (gridn(xval,yval,zval)*xp*yp*zp + gridn(xval,yval,zval+1)*xp*yp*zm + &
	gridn(xval+1,yval,zval)*xm*yp*zp + gridn(xval+1,yval,zval+1)*xm*yp*zm + &
	gridn(xval,yval+1,zval)*xp*ym*zp + gridn(xval,yval+1,zval+1)*xp*ym*zm + &
	gridn(xval+1,yval+1,zval)*xm*ym*zp + gridn(xval+1,yval+1,zval+1)*xm*ym*zm)

    inter_jx = (curl_x(xval,yval,zval)*xp*yp*zp + curl_x(xval,yval,zval+1)*xp*yp*zm + &
	curl_x(xval+1,yval,zval)*xm*yp*zp + curl_x(xval+1,yval,zval+1)*xm*yp*zm + &
	curl_x(xval,yval+1,zval)*xp*ym*zp + curl_x(xval,yval+1,zval+1)*xp*ym*zm + &
	curl_x(xval+1,yval+1,zval)*xm*ym*zp + curl_x(xval+1,yval+1,zval+1)*xm*ym*zm) 

    inter_jy = (curl_y(xval,yval,zval)*xp*yp*zp + curl_y(xval,yval,zval+1)*xp*yp*zm + &
	curl_y(xval+1,yval,zval)*xm*yp*zp + curl_y(xval+1,yval,zval+1)*xm*yp*zm + &
	curl_y(xval,yval+1,zval)*xp*ym*zp + curl_y(xval,yval+1,zval+1)*xp*ym*zm + &
	curl_y(xval+1,yval+1,zval)*xm*ym*zp + curl_y(xval+1,yval+1,zval+1)*xm*ym*zm)

    inter_jz = (curl_z(xval,yval,zval)*xp*yp*zp + curl_z(xval,yval,zval+1)*xp*yp*zm + &
	curl_z(xval+1,yval,zval)*xm*yp*zp + curl_z(xval+1,yval,zval+1)*xm*yp*zm + &
	curl_z(xval,yval+1,zval)*xp*ym*zp + curl_z(xval,yval+1,zval+1)*xp*ym*zm + &
	curl_z(xval+1,yval+1,zval)*xm*ym*zp + curl_z(xval+1,yval+1,zval+1)*xm*ym*zm)

    if (rand_num()*nsup .ge. real(inter_n,drk)) then
      np_i=np_i-1 
    else
      call maxwellian(real(T_i,drk), real(1.,drk), &
                      rvi(4,np_i), rvi(5,np_i), rvi(6,np_i))
#if (init_scheme == initharris)
call veloadd(rvi(4,np_i),rvi(5,np_i),rvi(6,np_i), T_i/(T_i+T_e)*&
	inter_jx/inter_n, T_i/(T_i+T_e)*inter_jy/inter_n, T_i/(T_i+T_e)*&
	inter_jz/inter_n )
#elif (init_scheme == initharris2)
      call veloadd(rvi(4,np_i),rvi(5,np_i),rvi(6,np_i), &
                   real( 2*T_i/w0*sign,drk), real(0.,drk), real(0.,drk) )
#endif
    endif
  enddo

! load electrons in current sheets
  do i=1,nx*ny*nz*ppg*nsup/n_0
    np_e=np_e+1
    if (np_e .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
    endif
! accept or reject particle according to local density
    call location(rve(1,np_e),rve(2,np_e),rve(3,np_e))

! See above comment
    rve(1,np_e) = nx*dx*(my_pex + rve(1,np_e)/lx)
    rve(2,np_e) = ny*dy*(my_pey + rve(2,np_e)/ly)
    rve(3,np_e) = nz*dz*(my_pez + rve(3,np_e)/lz)

! Interpolate density.  First find gridpoint _below_ point.
    xval = int(rve(1,np_e)/dx-my_pex*nx+0.5)
    yval = int(rve(2,np_e)/dy-my_pey*ny+0.5)
    zval = int(rve(3,np_e)/dz-my_pez*nz+0.5)
    xm=rve(1,np_e)/dx-my_pex*nx+0.5-xval; xp=1.-xm
    ym=rve(2,np_e)/dy-my_pey*ny+0.5-yval; yp=1.-ym
#if (nz == 1 && pez == 1)
    zm=1.; zp=0.
#else
    zm=rve(3,np_e)/dz-my_pez*nz+0.5-zval; zp=1.-zm
#endif
    inter_n = (gridn(xval,yval,zval)*xp*yp*zp + gridn(xval,yval,zval+1)*xp*yp*zm + &
	gridn(xval+1,yval,zval)*xm*yp*zp + gridn(xval+1,yval,zval+1)*xm*yp*zm + &
	gridn(xval,yval+1,zval)*xp*ym*zp + gridn(xval,yval+1,zval+1)*xp*ym*zm + &
	gridn(xval+1,yval+1,zval)*xm*ym*zp + gridn(xval+1,yval+1,zval+1)*xm*ym*zm)

    inter_jx = (curl_x(xval,yval,zval)*xp*yp*zp + curl_x(xval,yval,zval+1)*xp*yp*zm + &
	curl_x(xval+1,yval,zval)*xm*yp*zp + curl_x(xval+1,yval,zval+1)*xm*yp*zm + &
	curl_x(xval,yval+1,zval)*xp*ym*zp + curl_x(xval,yval+1,zval+1)*xp*ym*zm + &
	curl_x(xval+1,yval+1,zval)*xm*ym*zp + curl_x(xval+1,yval+1,zval+1)*xm*ym*zm) 

    inter_jy = (curl_y(xval,yval,zval)*xp*yp*zp + curl_y(xval,yval,zval+1)*xp*yp*zm + &
	curl_y(xval+1,yval,zval)*xm*yp*zp + curl_y(xval+1,yval,zval+1)*xm*yp*zm + &
	curl_y(xval,yval+1,zval)*xp*ym*zp + curl_y(xval,yval+1,zval+1)*xp*ym*zm + &
	curl_y(xval+1,yval+1,zval)*xm*ym*zp + curl_y(xval+1,yval+1,zval+1)*xm*ym*zm)

    inter_jz = (curl_z(xval,yval,zval)*xp*yp*zp + curl_z(xval,yval,zval+1)*xp*yp*zm + &
	curl_z(xval+1,yval,zval)*xm*yp*zp + curl_z(xval+1,yval,zval+1)*xm*yp*zm + &
	curl_z(xval,yval+1,zval)*xp*ym*zp + curl_z(xval,yval+1,zval+1)*xp*ym*zm + &
	curl_z(xval+1,yval+1,zval)*xm*ym*zp + curl_z(xval+1,yval+1,zval+1)*xm*ym*zm)

    if (rand_num()*nsup .ge. real(inter_n,drk)) then
      np_e=np_e-1 
    else
      call maxwellian(real(T_e,drk), real(m_e,drk), &
                      rve(4,np_e), rve(5,np_e), rve(6,np_e))
      sign=1; if (rve(2,np_e) .gt. ly/2.) sign=-1

! This kludge assumes the electrons and ions get the normal Harris
! sheet values, but the electrons get anything extra left over from, for
! instance, the perturbation.  Neater might have been to give the
! electrons all the current.

#if (init_scheme == initharris)
call veloadd(rve(4,np_e),rve(5,np_e),rve(6,np_e), -T_e/(T_i+T_e)*&
	inter_jx/inter_n,  -T_e/(T_i+T_e)*inter_jy/inter_n,  -T_e/(T_i+T_e)*&
	inter_jz/inter_n )

#elif (init_scheme == initharris2)
      call veloadd(rve(4,np_e),rve(5,np_e),rve(6,np_e), &
                     real(-2*T_e/w0*sign,drk), real(0.,drk), real(0.,drk) )
#endif
    endif
  enddo
  n_avg=ppg/n_0
#ifndef subtract_average_rho
! force global charge neutrality
  np_i=min(np_i,np_e); np_e=min(np_i,np_e)
#endif
  call redistribute(rvi,np_i)
  call redistribute(rve,np_e)

  call output(0,.true.)
End Subroutine init_harris_all
