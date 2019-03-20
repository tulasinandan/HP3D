!******************************************************************************
!             3D particle code: double harris sheet with perturbation
!                           Andreas Zeiler, 1999
!
!                      LATEST CHANGE: April 16, 2003
!
!******************************************************************************

! initperturb: current along z
! initperturb2: current along x

#include "param"
#define initperturb 1
#define initperturb2 2
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

Subroutine init_perturb_all()
  Use pe_env
  Use partfield
  Implicit None

  Real(kind=drk) gx,gy,gz,yy,b,n,rand_num,sech,nsup,j
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez),&
                          pi2=6.2831853, pi4=2.*pi2
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: curl_x,curl_y,curl_z,gridn
  Integer i,x,y,z,pe
  Real(kind=drk) :: realx,realy,realz
  Integer :: xval,yval,zval
  Real(kind=drk) :: xp,yp,zp,xm,ym,zm,inter_n,inter_jx,inter_jy,inter_jz

  gx(x)=((x-0.5)+my_pex*nx)*dx
  gy(y)=((y-0.5)+my_pey*ny)*dy
  gz(z)=((z-0.5)+my_pez*nz)*dz
  sech(yy)=1/cosh(yy)

! shape of the sheet
  b(yy)=tanh((yy-ly*.25)/w0)-tanh((yy-ly*.75)/w0)+&
	tanh((yy-ly*1.25)/w0)-tanh((yy+ly*.25)/w0)+1
  j(yy) = -(sech((yy-ly*.25)/w0)**2-sech((yy-ly*.75)/w0)**2+&
        sech((yy-ly*1.25)/w0)**2-sech((yy+ly*.25)/w0)**2)/w0
  n(yy) = (1-(tanh((yy-ly*.25)/w0)-tanh((yy-ly*.75)/w0)+&
      tanh((yy-ly*1.25)/w0)-tanh((yy+ly*.25)/w0)+1)**2)/(2*(T_i+T_e))
  nsup = 1.

! Anything processor dependent must be done after this next line.
  call init_pe_env()
  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** init: compiled for different number of PEs *****'
    call exitallpes()
  endif

! write parameters to log-file
  if (myproc==0) then
#if (init_scheme == initperturb)
    write(6,*) '********** perturb **********'
#elif (init_scheme == initperturb2)
    write(6,*) '********** initperturb2 *********'
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
#if (init_scheme == initperturb)
    write(6,*) '     psi0 = ',psi0
#endif
#ifdef relativistic
    write(6,*) '     relativistic initialization'
#endif
    write(6,*) '********** parameters **********'
  endif

! set electric and magnetic field (no kinetic equilibrium with guide field bz)
  do z=1,nz; do y=1,ny; do x=1,nx
#if (init_scheme == initperturb)
    b1x(x,y,z)=b(gy(y)); b1y(x,y,z)=0.; b1z(x,y,z)=b0
#elif (init_scheme == initperturb2)
    b1x(x,y,z)=-b0; b1y(x,y,z)=0.; b1z(x,y,z)=b(gy(y))
#endif
    e1x(x,y,z)=0.; e1y(x,y,z)=0.; e1z(x,y,z)=0.
  enddo; enddo; enddo 

! add perturbation to magnetic field to start reconnection
#ifdef bxpert
  do z=1,nz; do y=1,ny; do x=1,nx
    b1x(x,y,z)=b1x(x,y,z)-bxpert*sin(pi2*gx(x)/lx)*sin(pi4*gy(y)/ly)
    b1y(x,y,z)=b1y(x,y,z)+bxpert*ly/(2.*lx)*cos(pi2*gx(x)/lx)*(1-cos(pi4*gy(y)/ly))
  enddo; enddo; enddo 
#else
  do z=1,nz; do y=1,ny; do x=1,nx
    b1x(x,y,z)=b1x(x,y,z)-psi0*pi4/ly*sin(pi2*gx(x)/lx)*sin(pi4*gy(y)/ly)
    b1y(x,y,z)=b1y(x,y,z)+psi0*pi2/lx*cos(pi2*gx(x)/lx)*(1-cos(pi4*gy(y)/ly))
  enddo; enddo; enddo
#endif

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

  np_i=0; np_e=0
! load particles, background density.  Note that the perturbation is
! made in the background density.  If there is no such thing, J will not
! be curl B.
#ifndef no_background_density
  do i=1,nx*ny*nz*ppg
    np_i=np_i+1
    if (np_i .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
    endif

    call location(rvi(1,np_i),rvi(2,np_i),rvi(3,np_i))

! Move the particle to current processor.  Need to be careful
! that particles are not slightly shifted across boundaries due
! to arithmetic.  
    rvi(1,np_i) = lx*my_pex/pex + mod(rvi(1,np_i),lx/pex)
    rvi(2,np_i) = ly*my_pey/pey + mod(rvi(2,np_i),ly/pey)
    rvi(3,np_i) = lz*my_pez/pez + mod(rvi(3,np_i),lz/pez)
    if(rvi(1,np_i) .ge. lx*(my_pex+1)/pex) then
      rvi(1,np_i) = rvi(1,np_i) - lx/pex
    endif
    if(rvi(2,np_i) .ge. ly*(my_pey+1)/pey) then
      rvi(2,np_i) = rvi(2,np_i) - ly/pey
    endif
    if(rvi(3,np_i) .ge. lz*(my_pez+1)/pez) then
      rvi(3,np_i) = rvi(3,np_i) - lz/pez
    endif

! Interpolate density and current.  First find gridpoint _below_ point.
    xval = int(rvi(1,np_i)/dx-my_pex*nx+0.5)
    yval = int(rvi(2,np_i)/dy-my_pey*ny+0.5)
    zval = int(rvi(3,np_i)/dz-my_pez*nz+0.5)
    xm=rvi(1,np_i)/dx-my_pex*nx+0.5-xval; xp=1.-xm
    ym=rvi(2,np_i)/dy-my_pey*ny+0.5-yval; yp=1.-ym
    zm=rvi(3,np_i)/dz-my_pez*nz+0.5-zval; zp=1.-zm

    inter_jx = (curl_x(xval,yval,zval)*xp*yp*zp + &
      curl_x(xval,yval,zval+1)*xp*yp*zm + &
      curl_x(xval+1,yval,zval)*xm*yp*zp + &
      curl_x(xval+1,yval,zval+1)*xm*yp*zm + &
      curl_x(xval,yval+1,zval)*xp*ym*zp + &
      curl_x(xval,yval+1,zval+1)*xp*ym*zm + &
      curl_x(xval+1,yval+1,zval)*xm*ym*zp + &
      curl_x(xval+1,yval+1,zval+1)*xm*ym*zm) 

    inter_jy = (curl_y(xval,yval,zval)*xp*yp*zp + &
      curl_y(xval,yval,zval+1)*xp*yp*zm + &
      curl_y(xval+1,yval,zval)*xm*yp*zp + &
      curl_y(xval+1,yval,zval+1)*xm*yp*zm + &
      curl_y(xval,yval+1,zval)*xp*ym*zp + &
      curl_y(xval,yval+1,zval+1)*xp*ym*zm + &
      curl_y(xval+1,yval+1,zval)*xm*ym*zp + &
      curl_y(xval+1,yval+1,zval+1)*xm*ym*zm)

    inter_jz = (curl_z(xval,yval,zval)*xp*yp*zp + &
      curl_z(xval,yval,zval+1)*xp*yp*zm + &
      curl_z(xval+1,yval,zval)*xm*yp*zp + & 
      curl_z(xval+1,yval,zval+1)*xm*yp*zm + &
      curl_z(xval,yval+1,zval)*xp*ym*zp + & 
      curl_z(xval,yval+1,zval+1)*xp*ym*zm + &
      curl_z(xval+1,yval+1,zval)*xm*ym*zp + & 
      curl_z(xval+1,yval+1,zval+1)*xm*ym*zm)

    call maxwellian(real(T_i,drk), real(1.,drk), &
                    rvi(4,np_i), rvi(5,np_i), rvi(6,np_i))
    call veloadd(rvi(4,np_i),rvi(5,np_i),rvi(6,np_i), T_i/(T_i+T_e)*&
	inter_jx/n_0, T_i/(T_i+T_e)*inter_jy/n_0, &
	T_i/(T_i+T_e)*(inter_jz-j(real(rvi(2,np_i),drk)))/n_0 )
  enddo

  do i=1,nx*ny*nz*ppg
    np_e=np_e+1
    if (np_e .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
    endif
    call location(rve(1,np_e),rve(2,np_e),rve(3,np_e))
! Move the particle to current processor.  Need to be careful
! that particles are not slightly shifted across boundaries due
! to arithmetic. 
    rve(1,np_e) = lx*my_pex/pex + mod(rve(1,np_e),lx/pex)
    rve(2,np_e) = ly*my_pey/pey + mod(rve(2,np_e),ly/pey)
    rve(3,np_e) = lz*my_pez/pez + mod(rve(3,np_e),lz/pez)
    if(rve(1,np_e) .ge. lx*(my_pex+1)/pex) then
      rve(1,np_e) = rve(1,np_e) - lx/pex
    endif
    if(rve(2,np_e) .ge. ly*(my_pey+1)/pey) then
      rve(2,np_e) = rve(2,np_e) - ly/pey
    endif
    if(rve(3,np_e) .ge. lz*(my_pez+1)/pez) then
      rve(3,np_e) = rve(3,np_e) - lz/pez
    endif

! Interpolate density.  First find gridpoint _below_ point.
    xval = int(rve(1,np_e)/dx-my_pex*nx+0.5)
    yval = int(rve(2,np_e)/dy-my_pey*ny+0.5)
    zval = int(rve(3,np_e)/dz-my_pez*nz+0.5)
    xm=rve(1,np_e)/dx-my_pex*nx+0.5-xval; xp=1.-xm
    ym=rve(2,np_e)/dy-my_pey*ny+0.5-yval; yp=1.-ym
    zm=rve(3,np_e)/dz-my_pez*nz+0.5-zval; zp=1.-zm

    inter_jx = (curl_x(xval,yval,zval)*xp*yp*zp + &
      curl_x(xval,yval,zval+1)*xp*yp*zm + &
      curl_x(xval+1,yval,zval)*xm*yp*zp + &
      curl_x(xval+1,yval,zval+1)*xm*yp*zm + &
      curl_x(xval,yval+1,zval)*xp*ym*zp + &
      curl_x(xval,yval+1,zval+1)*xp*ym*zm + &
      curl_x(xval+1,yval+1,zval)*xm*ym*zp + &
      curl_x(xval+1,yval+1,zval+1)*xm*ym*zm) 

    inter_jy = (curl_y(xval,yval,zval)*xp*yp*zp + &
      curl_y(xval,yval,zval+1)*xp*yp*zm + &
      curl_y(xval+1,yval,zval)*xm*yp*zp + &
      curl_y(xval+1,yval,zval+1)*xm*yp*zm + &
      curl_y(xval,yval+1,zval)*xp*ym*zp + &
      curl_y(xval,yval+1,zval+1)*xp*ym*zm + &
      curl_y(xval+1,yval+1,zval)*xm*ym*zp + &
      curl_y(xval+1,yval+1,zval+1)*xm*ym*zm)

    inter_jz = (curl_z(xval,yval,zval)*xp*yp*zp + &
      curl_z(xval,yval,zval+1)*xp*yp*zm + &
      curl_z(xval+1,yval,zval)*xm*yp*zp + &
      curl_z(xval+1,yval,zval+1)*xm*yp*zm + &
      curl_z(xval,yval+1,zval)*xp*ym*zp + &
      curl_z(xval,yval+1,zval+1)*xp*ym*zm + &
      curl_z(xval+1,yval+1,zval)*xm*ym*zp + &
      curl_z(xval+1,yval+1,zval+1)*xm*ym*zm)

    call maxwellian(real(T_e,drk), real(m_e,drk), &
                    rve(4,np_e), rve(5,np_e), rve(6,np_e))
    call veloadd(rve(4,np_e),rve(5,np_e),rve(6,np_e), -T_e/(T_i+T_e)*&
      inter_jx/n_0,  -T_e/(T_i+T_e)*inter_jy/n_0,  -T_e/(T_i+T_e)*&
      (inter_jz-j(real(rve(2,np_e),drk)))/n_0 )
  enddo
#endif

! load ions in current sheets
  do i=1,nx*ny*nz*ppg*nsup/n_0
    np_i=np_i+1
    if (np_i .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
    endif
! accept or reject particle according to local density
    call location(rvi(1,np_i),rvi(2,np_i),rvi(3,np_i))

! Move the particle to current processor.  Need to be careful
! that particles are not slightly shifted across boundaries due
! to arithmetic.  
    rvi(1,np_i) = lx*my_pex/pex + mod(rvi(1,np_i),lx/pex)
    rvi(2,np_i) = ly*my_pey/pey + mod(rvi(2,np_i),ly/pey)
    rvi(3,np_i) = lz*my_pez/pez + mod(rvi(3,np_i),lz/pez)
    if(rvi(1,np_i) .ge. lx*(my_pex+1)/pex) then
      rvi(1,np_i) = rvi(1,np_i) - lx/pex
    endif
    if(rvi(2,np_i) .ge. ly*(my_pey+1)/pey) then
      rvi(2,np_i) = rvi(2,np_i) - ly/pey
    endif
    if(rvi(3,np_i) .ge. lz*(my_pez+1)/pez) then
      rvi(3,np_i) = rvi(3,np_i) - lz/pez
    endif

    if (rand_num()*nsup .ge. n(real(rvi(2,np_i),drk))) then
      np_i=np_i-1 
    else
      call maxwellian(real(T_i,drk), real(1.,drk), &
                      rvi(4,np_i), rvi(5,np_i), rvi(6,np_i))
#if (init_scheme == initperturb)
      call veloadd(rvi(4,np_i),rvi(5,np_i),rvi(6,np_i), &
        real(0.,drk),real(0.,drk), T_i/(T_i+T_e)* &
        j(real(rvi(2,np_i),drk))/n(real(rvi(2,np_i),drk)))
#elif (init_scheme == initperturb2)
      call veloadd(rvi(4,np_i),rvi(5,np_i),rvi(6,np_i), &
        T_i/(T_i+T_e)*j(real(rvi(2,np_i),drk))/n(real(rvi(2,np_i),drk)), &
        real(0.,drk), real(0.,drk) )
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

! Move the particle to current processor.  Need to be careful
! that particles are not slightly shifted across boundaries due
! to arithmetic. 
    rve(1,np_e) = lx*my_pex/pex + mod(rve(1,np_e),lx/pex)
    rve(2,np_e) = ly*my_pey/pey + mod(rve(2,np_e),ly/pey)
    rve(3,np_e) = lz*my_pez/pez + mod(rve(3,np_e),lz/pez)
    if(rve(1,np_e) .ge. lx*(my_pex+1)/pex) then
      rve(1,np_e) = rve(1,np_e) - lx/pex
    endif
    if(rve(2,np_e) .ge. ly*(my_pey+1)/pey) then
      rve(2,np_e) = rve(2,np_e) - ly/pey
    endif
    if(rve(3,np_e) .ge. lz*(my_pez+1)/pez) then
      rve(3,np_e) = rve(3,np_e) - lz/pez
    endif

    if (rand_num()*nsup .ge. n(real(rve(2,np_e),drk))) then
      np_e=np_e-1 
    else
      call maxwellian(real(T_e,drk), real(m_e,drk), &
                      rve(4,np_e), rve(5,np_e), rve(6,np_e))
#if (init_scheme == initperturb)
      call veloadd(rve(4,np_e),rve(5,np_e),rve(6,np_e), &
        real(0.,drk), real(0.,drk), -T_e/(T_i+T_e)*&
        j(real(rve(2,np_e),drk))/n(real(rve(2,np_e),drk)))
#elif (init_scheme == initperturb2)
      call veloadd(rve(4,np_e),rve(5,np_e),rve(6,np_e), &
        -T_e/(T_i+T_e)*j(real(rve(2,np_e),drk))/n(real(rve(2,np_e),drk)), &
        real(0.,drk), real(0.,drk) )
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
End Subroutine init_perturb_all
