!*****************************************************************************
!                   3D particle code: particle dynamics
!                           Andreas Zeiler, 1999
!
!                        LATEST CHANGE: May 3, 2001
!
!*****************************************************************************

#include "param"

!-----------------------------------------------------------------------------
!                      step particles: positions
!-----------------------------------------------------------------------------

Subroutine stepx(rv,np)
! steps particle species in space according to velocity v by dt/2 
! periodic boundary conditions in x,z; variable in y
! rv vector with positions and velocities, np number of particles on PE

  Use pe_env
  Use partfield
#ifdef EXPANSION
  Use expansion_params
#endif
  Use timer
  Implicit none
  Real(kind=prk), Dimension(6,maxparticles) :: rv
!#ifdef EXPANSION
!  Real(kind=prk), Dimension(3,maxparticles) :: velocity
!#endif
  Real(kind=prk), Parameter :: epsilon=1e-4
  Integer :: np,i,ix,iy,iz
  Logical :: non_periodic

  call start_timer()
!#ifdef EXPANSION
!   if (myproc == 0) write(6,*) 'ww, dwow, aa, daoa', ww, dwow, aa, daoa
!   velocity(1,:)=rv(4,:)/ww
!   velocity(2,:)=rv(5,:)/aa
!   velocity(3,:)=rv(6,:)/aa
!#endif

  if (non_periodic()) then

! MS: Moved within the non_periodic loop
! check which components of momentum need to be reverted at boundary
  call bound_particle(ix,iy,iz)

    do i=1, np
#ifdef EXPANSION
!! NEED TO ADD THE MIRROR FORCE STILL
!     rv(1:3,i)=rv(1:3,i)+(dt)/2.*velocity(1:3,i)
      rv(1,i)=rv(1,i)+dt*rv(4,i)/(2.*ww)
      rv(2,i)=rv(2,i)+dt*rv(5,i)/(2.*aa)
      rv(3,i)=rv(3,i)+dt*rv(6,i)/(2.*aa)
#else
      rv(1:3,i)=rv(1:3,i)+(dt)/2.*rv(4:6,i)
#endif
      do while(rv(2,i) .lt. epsilon*ly .or. rv(2,i) .gt. (1-epsilon)*ly)
        if (rv(2,i) .le. epsilon*ly) then
          rv(2,i)=2*epsilon*ly-rv(2,i)
        else
          rv(2,i)=2*(1-epsilon)*ly-rv(2,i)
        endif
        rv(4,i)=ix*rv(4,i); rv(5,i)=iy*rv(5,i); rv(6,i)=iz*rv(6,i)
      enddo
      if(rv(1,i)<0. .or. rv(1,i)>=lx) rv(1,i)=mod(rv(1,i)+lx,real(lx,prk))
      if(rv(3,i)<0. .or. rv(3,i)>=lz) rv(3,i)=mod(rv(3,i)+lz,real(lz,prk))
    enddo
  else

! Removed a +100*lx,ly,lz from the position step and added if's around
! the mod step.  The former was to make sure position was never
! negative; that's now done in the mod step.  Note that the F90 version
! of mod behaves weirdly for negative numbers.  Rather than the somewhat
! inelegant +lx I could use modulo, but the latter is much slower and 
! prone to round-off problems.

! M. Shay, 8/20/2004
! Modified the particle movement below for the case of a quasi-1D
     ! system where ny = 2. In this system, it is possible for the a
     ! particle to travel more than ly away from the edge of a box in
     ! one particle time step. Done the old way, this particle would
     ! not be placed in the active grids when it went more than 1 ly
     ! in the negative direction. The current way with two mods fixes
     ! that. 

! M. Shay, 5/15/2005
! As mentioned above, when the code is run in quasi-1D mode there can
! problems with particles not getting put onto the active grid. If
! the code is run in quasi-2D mode, the modified function with two
! mod's is used instead. Quasi-2D means a total of 2 grid points
! along the y direction across the WHOLE SYSTEM (not just across one
! processor). 

    do i=1, np
#ifdef EXPANSION
!! NEED TO ADD THE MIRROR FORCE STILL
!     rv(1:3,i)=rv(1:3,i)+(dt)/2.*velocity(1:3,i)
      rv(1,i)=rv(1,i)+dt*rv(4,i)/(2.*ww)
      rv(2,i)=rv(2,i)+dt*rv(5,i)/(2.*aa)
      rv(3,i)=rv(3,i)+dt*rv(6,i)/(2.*aa)
#else
      rv(1:3,i)=rv(1:3,i)+(dt)/2.*rv(4:6,i)
#endif
      if(rv(1,i)<0. .or. rv(1,i)>=lx) rv(1,i)=mod(rv(1,i)+lx,real(lx,prk))
#     if (ny==2 && pey==1)
        if(rv(2,i)<0. .or. rv(2,i)>=ly) rv(2,i)=mod(mod(rv(2,i),real(ly,prk))+ly,real(ly,prk))
#     else
        if(rv(2,i)<0. .or. rv(2,i)>=ly) rv(2,i)=mod(rv(2,i)+ly,real(ly,prk))
#     endif
      if(rv(3,i)<0. .or. rv(3,i)>=lz) rv(3,i)=mod(rv(3,i)+lz,real(lz,prk))
    enddo

  endif
!  call stop_timer(t_stepx,n_stepx)
! transfer particles, that cross the boundary, to appropriate PE
  call redistribute(rv,np)

  call stop_timer(t_stepx,n_stepx)

End Subroutine stepx

!-----------------------------------------------------------------------------
!                      step particles: velocities
!-----------------------------------------------------------------------------

Subroutine stepv(rv,np,mass)
! steps particle velocities by dt 
! rv vector with positions and velocities, np number of particles on PE
! mass particle mass (may be negative to account for sign of charge)

  Use pe_env
  Use partfield
#ifdef EXPANSION
  Use expansion_params
#endif
  Use timer
  Implicit none
  Real(kind=prk), Dimension(6,maxparticles) :: rv
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez)
  Integer :: np,i,x,y,z,m,n
  Real(kind=drk) :: mass,ex,ey,ez,bx,by,bz,xp,yp,zp,xm,ym,zm,a,gamma
  Real(kind=drk), Dimension(3,3) :: a1,a2
  Real(kind=drk), Dimension(3) :: tmp

  call start_timer()
! smooth electric field, update boundary layers
  call smooth_e(e1x,e1y,e1z,esx,esy,esz)
  call bound_e(esx,esy,esz)
  call bound_b(b1x,b1y,b1z) 

! step particle velocity
#if (nz == 1 && pez == 1) 
  z = 0; zm=1.; zp=0.
  do i=1, np
! interpolate electric and magnetic field
    x=int(rv(1,i)/dx-my_pex*nx+0.5)
    y=int(rv(2,i)/dy-my_pey*ny+0.5)
    xm=rv(1,i)/dx-my_pex*nx+0.5-x; xp=1.-xm
    ym=rv(2,i)/dy-my_pey*ny+0.5-y; yp=1.-ym
#else
  do i=1, np
! interpolate electric and magnetic field
    x=int(rv(1,i)/dx-my_pex*nx+0.5)
    y=int(rv(2,i)/dy-my_pey*ny+0.5)
    z=int(rv(3,i)/dz-my_pez*nz+0.5)
    xm=rv(1,i)/dx-my_pex*nx+0.5-x; xp=1.-xm
    ym=rv(2,i)/dy-my_pey*ny+0.5-y; yp=1.-ym
    zm=rv(3,i)/dz-my_pez*nz+0.5-z; zp=1.-zm
#endif

    ex = ( esx(x  ,y  ,z  )*xp*yp*zp + esx(x  ,y  ,z+1)*xp*yp*zm + &
           esx(x+1,y  ,z  )*xm*yp*zp + esx(x+1,y  ,z+1)*xm*yp*zm + &
           esx(x  ,y+1,z  )*xp*ym*zp + esx(x  ,y+1,z+1)*xp*ym*zm + &
           esx(x+1,y+1,z  )*xm*ym*zp + esx(x+1,y+1,z+1)*xm*ym*zm )*(dt)/(2*mass)
    ey = ( esy(x  ,y  ,z  )*xp*yp*zp + esy(x  ,y  ,z+1)*xp*yp*zm + &
           esy(x+1,y  ,z  )*xm*yp*zp + esy(x+1,y  ,z+1)*xm*yp*zm + &
           esy(x  ,y+1,z  )*xp*ym*zp + esy(x  ,y+1,z+1)*xp*ym*zm + &
           esy(x+1,y+1,z  )*xm*ym*zp + esy(x+1,y+1,z+1)*xm*ym*zm )*(dt)/(2*mass)
    ez = ( esz(x  ,y  ,z  )*xp*yp*zp + esz(x  ,y  ,z+1)*xp*yp*zm + &
           esz(x+1,y  ,z  )*xm*yp*zp + esz(x+1,y  ,z+1)*xm*yp*zm + &
           esz(x  ,y+1,z  )*xp*ym*zp + esz(x  ,y+1,z+1)*xp*ym*zm + &
           esz(x+1,y+1,z  )*xm*ym*zp + esz(x+1,y+1,z+1)*xm*ym*zm )*(dt)/(2*mass)
    bx = ( b1x(x  ,y  ,z  )*xp*yp*zp + b1x(x  ,y  ,z+1)*xp*yp*zm + &
           b1x(x+1,y  ,z  )*xm*yp*zp + b1x(x+1,y  ,z+1)*xm*yp*zm + &
           b1x(x  ,y+1,z  )*xp*ym*zp + b1x(x  ,y+1,z+1)*xp*ym*zm + &
           b1x(x+1,y+1,z  )*xm*ym*zp + b1x(x+1,y+1,z+1)*xm*ym*zm )*(dt)/(2*mass)
    by = ( b1y(x  ,y  ,z  )*xp*yp*zp + b1y(x  ,y  ,z+1)*xp*yp*zm + &
           b1y(x+1,y  ,z  )*xm*yp*zp + b1y(x+1,y  ,z+1)*xm*yp*zm + &
           b1y(x  ,y+1,z  )*xp*ym*zp + b1y(x  ,y+1,z+1)*xp*ym*zm + &
           b1y(x+1,y+1,z  )*xm*ym*zp + b1y(x+1,y+1,z+1)*xm*ym*zm )*(dt)/(2*mass)
    bz = ( b1z(x  ,y  ,z  )*xp*yp*zp + b1z(x  ,y  ,z+1)*xp*yp*zm + &
           b1z(x+1,y  ,z  )*xm*yp*zp + b1z(x+1,y  ,z+1)*xm*yp*zm + &
           b1z(x  ,y+1,z  )*xp*ym*zp + b1z(x  ,y+1,z+1)*xp*ym*zm + &
           b1z(x+1,y+1,z  )*xm*ym*zp + b1z(x+1,y+1,z+1)*xm*ym*zm )*(dt)/(2*mass)
! first half of electric impulse
! relativistic case: convert velocity into relativistic momentum
#ifdef relativistic
    gamma=1/sqrt(1-(rv(4,i)**2+rv(5,i)**2+rv(6,i)**2)/c_2)
#else
    gamma=1.
#endif

#ifdef EXPANSION
    rv(4,i) = rv(4,i)*gamma+ex-dwow*rv(4,i)*dt/(2*mass)
    rv(5,i) = rv(5,i)*gamma+ey-daoa*rv(5,i)*dt/(2*mass)
    rv(6,i) = rv(6,i)*gamma+ez-daoa*rv(6,i)*dt/(2*mass)
#else
    rv(4,i) = rv(4,i)*gamma+ex
    rv(5,i) = rv(5,i)*gamma+ey
    rv(6,i) = rv(6,i)*gamma+ez
#endif
! magnetic rotation, matrix A: transpose a1, inverse a2/a
#ifdef relativistic
! reduce rotation angle by factor gamma
    gamma=sqrt(1+(rv(4,i)**2+rv(5,i)**2+rv(6,i)**2)/c_2)
    bx=bx/gamma
    by=by/gamma
    bz=bz/gamma
#endif
!   a1(1,1)=  1; a1(1,2)= bz; a1(1,3) =-by
!   a1(2,1)=-bz; a1(2,2)=  1; a1(2,3) = bx
!   a1(3,1)= by; a1(3,2)=-bx; a1(3,3) =  1
!   a2(1,1)= 1+bx**2; a2(1,2)=bx*by+bz; a2(1,3)=bx*bz-by
!   a2(2,1)=bx*by-bz; a2(2,2)= 1+by**2; a2(2,3)=by*bz+bx
!   a2(3,1)=bx*bz+by; a2(3,2)=by*bz-bx; a2(3,3)= 1+bz**2
    a=1+bx**2+by**2+bz**2
! rotate first with transpose, then with inverse
!    do m=1,3; tmp(m)=0.; enddo
!    do m=1,3; do n=1,3
!      tmp(m)=tmp(m)+a1(m,n)*rv(n+3,i)
!    enddo; enddo  
!    do m=1,3; rv(m+3,i)=0.; enddo
!    do m=1,3; do n=1,3
!      rv(m+3,i)=rv(m+3,i)+a2(m,n)/a*tmp(n)
!   enddo; enddo  
!tmp(1) = rv(4,i)*(1+bx**2-by**2-bz**2)+rv(5,i)*2*(by*bx+bz)+rv(6,i)*2*(bx*bz-by)
!tmp(2) = rv(4,i)*2*(by*bx-bz)+rv(5,i)*(1-bx**2+by**2-bz**2)+rv(6,i)*2*(by*bz+bx)
!tmp(3) = rv(4,i)*2*(bx*bz+by)+rv(5,i)*2*(by*bz-bx)+rv(6,i)*(1-bx**2-by**2+bz**2)
!rv(4:6,i) = tmp/a
tmp(1) = rv(4,i)+rv(5,i)*bz-rv(6,i)*by
tmp(2) = rv(5,i)+rv(6,i)*bx-rv(4,i)*bz
tmp(3) = rv(6,i)+rv(4,i)*by-rv(5,i)*bx
tmp(:)=tmp(:)*2/a
rv(4,i) = rv(4,i)+tmp(2)*bz-tmp(3)*by
rv(5,i) = rv(5,i)+tmp(3)*bx-tmp(1)*bz
rv(6,i) = rv(6,i)+tmp(1)*by-tmp(2)*bx

! second half of electric impulse
#ifdef EXPANSION
    rv(4,i) = rv(4,i)+ex-dwow*rv(4,i)*dt/(2*mass)
    rv(5,i) = rv(5,i)+ey-daoa*rv(5,i)*dt/(2*mass)
    rv(6,i) = rv(6,i)+ez-daoa*rv(6,i)*dt/(2*mass)
#else
    rv(4,i) = rv(4,i)+ex
    rv(5,i) = rv(5,i)+ey
    rv(6,i) = rv(6,i)+ez
#endif
#ifdef relativistic
! convert relativistic momentum into velocity
    gamma=sqrt(1+(rv(4,i)**2+rv(5,i)**2+rv(6,i)**2)/c_2)
    rv(4,i) = rv(4,i)/gamma
    rv(5,i) = rv(5,i)/gamma
    rv(6,i) = rv(6,i)/gamma
#endif
  enddo
  call stop_timer(t_stepv,n_stepv)
End Subroutine stepv

!-----------------------------------------------------------------------------
!                      calculate charge density
!-----------------------------------------------------------------------------

Subroutine calc_rho(rv,np,charge,reset)
! calculates charge density due to particle species with charge 'charge'
! rv vector with positions and velocities, np number of particles on PE
! charge is added to rho(x,y,z), rho is initialized to zero if reset=.true.

  Use pe_env
  Use partfield
  Use timer
#ifdef EXPANSION
  Use expansion_params
#endif
  Implicit none
  Real(kind=prk), Dimension(6,maxparticles) :: rv
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez)
  Real(kind=drk) :: charge, xm, ym, zm, xp, yp, zp, tmp
  Integer :: np,i,x,y,z,mc
  Logical :: reset

  call start_timer()
! reset array for charge density
  if (reset) then
    do x=0,nx+1; do y=0,ny+1; do z=0,nz+1
      rho(x,y,z)=0.
    enddo; enddo; enddo
!rho = 0.
  else
! reset boundary layers only
    do x=0,nx+1; do y=0,ny+1
      rho(x,y,0)=0.; rho(x,y,nz+1)=0.
    enddo; enddo
!rho(0:nx+1,0:ny+1,0)=0.; rho(0:nx+1,0:ny+1,nz+1)=0.
    do x=0,nx+1; do z=0,nz+1
      rho(x,0,z)=0.; rho(x,ny+1,z)=0.
    enddo; enddo
!rho(0:nx+1,0,0:nz+1)=0.; rho(0:nx+1,ny+1,0:nz+1)=0.
    do y=0,ny+1; do z=0,nz+1
      rho(0,y,z)=0.; rho(nx+1,y,z)=0.
    enddo; enddo
!rho(0,0:ny+1,0:nz+1)=0.; rho(nx+1,0:ny+1,0:nz+1)=0.
  endif

  call virtual_charge(mc)
! calculate charge density
#if (nz == 1 && pez == 1) 
  z = 0; zm=1.; zp=0.
  do i=1, np
    x=int(rv(1,i)/dx-my_pex*nx+0.5)
    y=int(rv(2,i)/dy-my_pey*ny+0.5)
    xm=rv(1,i)/dx-my_pex*nx+0.5-x; xp=1.-xm
    ym=rv(2,i)/dy-my_pey*ny+0.5-y; yp=1.-ym
#else
  do i=1, np
    x=int(rv(1,i)/dx-my_pex*nx+0.5)
    y=int(rv(2,i)/dy-my_pey*ny+0.5)
    z=int(rv(3,i)/dz-my_pez*nz+0.5)
    xm=rv(1,i)/dx-my_pex*nx+0.5-x; xp=1.-xm
    ym=rv(2,i)/dy-my_pey*ny+0.5-y; yp=1.-ym
    zm=rv(3,i)/dz-my_pez*nz+0.5-z; zp=1.-zm
#endif
    rho(x  ,y  ,z  )=rho(x  ,y  ,z  )+charge*xp*yp*zp/n_avg
    rho(x+1,y  ,z  )=rho(x+1,y  ,z  )+charge*xm*yp*zp/n_avg
    rho(x  ,y+1,z  )=rho(x  ,y+1,z  )+charge*xp*ym*zp/n_avg
    rho(x  ,y  ,z+1)=rho(x  ,y  ,z+1)+charge*xp*yp*zm/n_avg
    rho(x  ,y+1,z+1)=rho(x  ,y+1,z+1)+charge*xp*ym*zm/n_avg
    rho(x+1,y  ,z+1)=rho(x+1,y  ,z+1)+charge*xm*yp*zm/n_avg
    rho(x+1,y+1,z  )=rho(x+1,y+1,z  )+charge*xm*ym*zp/n_avg
    rho(x+1,y+1,z+1)=rho(x+1,y+1,z+1)+charge*xm*ym*zm/n_avg
    if (mc .ne. 0) then
      if ((my_pey==0 .and. y==0) .or. (my_pey==pey-1 .and. y==ny)) then
        tmp=yp; yp=ym; ym=tmp
        rho(x  ,y  ,z  )=rho(x  ,y  ,z  )+charge*xp*yp*zp/n_avg*mc
        rho(x+1,y  ,z  )=rho(x+1,y  ,z  )+charge*xm*yp*zp/n_avg*mc
        rho(x  ,y+1,z  )=rho(x  ,y+1,z  )+charge*xp*ym*zp/n_avg*mc
        rho(x  ,y  ,z+1)=rho(x  ,y  ,z+1)+charge*xp*yp*zm/n_avg*mc
        rho(x  ,y+1,z+1)=rho(x  ,y+1,z+1)+charge*xp*ym*zm/n_avg*mc
        rho(x+1,y  ,z+1)=rho(x+1,y  ,z+1)+charge*xm*yp*zm/n_avg*mc
        rho(x+1,y+1,z  )=rho(x+1,y+1,z  )+charge*xm*ym*zp/n_avg*mc
        rho(x+1,y+1,z+1)=rho(x+1,y+1,z+1)+charge*xm*ym*zm/n_avg*mc
      endif
    endif
  enddo
#ifdef EXPANSION
  rho = rho/aa**2
#endif
! transfer charge in boundary layer to adjacent PEs
  call bound_rho_j(rho)
  call stop_timer(t_calc_rho,n_calc_rho)
End Subroutine calc_rho

!-----------------------------------------------------------------------------
!                      calculate current density
!-----------------------------------------------------------------------------

Subroutine calc_j(rv,np,charge,reset)
! calculates current density due to particle species with charge 'charge'
! rv vector with positions and velocities, np number of particles on PE
! current is added to jx,y,z(x,y,z), initialized to zero if reset=.true.

  Use pe_env
  Use partfield
  Use timer
  Implicit none
  Real(kind=prk), Dimension(6,maxparticles) :: rv
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez)
  Real(kind=drk) :: charge, xm, ym, zm, xp, yp, zp, tmp
  Integer :: np,i,x,y,z,mcx,mcy,mcz
  Logical :: reset
  
  call start_timer()
! reset arrays for current density
  if (reset) then
    do x=0,nx+1; do y=0,ny+1; do z=0,nz+1
      jx(x,y,z)=0.; jy(x,y,z)=0.; jz(x,y,z)=0.
    enddo; enddo; enddo
!jx =0.; jy = 0.; jz = 0.
  else
! reset boundary layers only
    do x=0,nx+1; do y=0,ny+1
      jx(x,y,0)=0.; jx(x,y,nz+1)=0.
      jy(x,y,0)=0.; jy(x,y,nz+1)=0.
      jz(x,y,0)=0.; jz(x,y,nz+1)=0.
    enddo; enddo
!jx(0:nx+1,0:ny+1,0)=0.; jy(0:nx+1,0:ny+1,0)=0.; jz(0:nx+1,0:ny+1,0)=0.
!jx(0:nx+1,0:ny+1,nz+1)=0.; jy(0:nx+1,0:ny+1,nz+1)=0.; jz(0:nx+1,0:ny+1,nz+1)=0.
    do x=0,nx+1; do z=0,nz+1
      jx(x,0,z)=0.; jx(x,ny+1,z)=0.
      jy(x,0,z)=0.; jy(x,ny+1,z)=0.
      jz(x,0,z)=0.; jz(x,ny+1,z)=0.
    enddo; enddo
!jx(0:nx+1,0,0:nz+1)=0.; jy(0:nx+1,0,0:nz+1)=0.; jz(0:nx+1,0,0:nz+1)=0.
!jx(0:nx+1,ny+1,0:nz+1)=0.; jy(0:nx+1,ny+1,0:nz+1)=0.; jz(0:nx+1,ny+1,0:nz+1)=0.
    do y=0,ny+1; do z=0,nz+1
      jx(0,y,z)=0.; jx(nx+1,y,z)=0.
      jy(0,y,z)=0.; jy(nx+1,y,z)=0.
      jz(0,y,z)=0.; jz(nx+1,y,z)=0.
    enddo; enddo
!jx(0,0:ny+1,0:nz+1)=0.; jy(0,0:ny+1,0:nz+1)=0.; jz(0,0:ny+1,0:nz+1)=0.
!jx(nx+1,0:ny+1,0:nz+1)=0.; jy(nx+1,0:ny+1,0:nz+1)=0.; jz(nx+1,0:ny+1,0:nz+1)=0.
  endif

  call virtual_current(mcx,mcy,mcz)
! calculate current density
#if (nz == 1 && pez == 1) 
  z = 0; zm=1.; zp=0.
  do i=1, np
    x=int(rv(1,i)/dx-my_pex*nx+0.5)
    y=int(rv(2,i)/dy-my_pey*ny+0.5)
    xm=rv(1,i)/dx-my_pex*nx+0.5-x; xp=1.-xm
    ym=rv(2,i)/dy-my_pey*ny+0.5-y; yp=1.-ym
#else
  do i=1, np
    x=int(rv(1,i)/dx-my_pex*nx+0.5)
    y=int(rv(2,i)/dy-my_pey*ny+0.5)
    z=int(rv(3,i)/dz-my_pez*nz+0.5)
    xm=rv(1,i)/dx-my_pex*nx+0.5-x; xp=1.-xm
    ym=rv(2,i)/dy-my_pey*ny+0.5-y; yp=1.-ym
    zm=rv(3,i)/dz-my_pez*nz+0.5-z; zp=1.-zm
#endif
    jx(x  ,y  ,z  )=jx(x  ,y  ,z  )+charge*xp*yp*zp*rv(4,i)/n_avg
    jx(x+1,y  ,z  )=jx(x+1,y  ,z  )+charge*xm*yp*zp*rv(4,i)/n_avg
    jx(x  ,y+1,z  )=jx(x  ,y+1,z  )+charge*xp*ym*zp*rv(4,i)/n_avg
    jx(x  ,y  ,z+1)=jx(x  ,y  ,z+1)+charge*xp*yp*zm*rv(4,i)/n_avg
    jx(x  ,y+1,z+1)=jx(x  ,y+1,z+1)+charge*xp*ym*zm*rv(4,i)/n_avg
    jx(x+1,y  ,z+1)=jx(x+1,y  ,z+1)+charge*xm*yp*zm*rv(4,i)/n_avg
    jx(x+1,y+1,z  )=jx(x+1,y+1,z  )+charge*xm*ym*zp*rv(4,i)/n_avg
    jx(x+1,y+1,z+1)=jx(x+1,y+1,z+1)+charge*xm*ym*zm*rv(4,i)/n_avg
    jy(x  ,y  ,z  )=jy(x  ,y  ,z  )+charge*xp*yp*zp*rv(5,i)/n_avg
    jy(x+1,y  ,z  )=jy(x+1,y  ,z  )+charge*xm*yp*zp*rv(5,i)/n_avg
    jy(x  ,y+1,z  )=jy(x  ,y+1,z  )+charge*xp*ym*zp*rv(5,i)/n_avg
    jy(x  ,y  ,z+1)=jy(x  ,y  ,z+1)+charge*xp*yp*zm*rv(5,i)/n_avg
    jy(x  ,y+1,z+1)=jy(x  ,y+1,z+1)+charge*xp*ym*zm*rv(5,i)/n_avg
    jy(x+1,y  ,z+1)=jy(x+1,y  ,z+1)+charge*xm*yp*zm*rv(5,i)/n_avg
    jy(x+1,y+1,z  )=jy(x+1,y+1,z  )+charge*xm*ym*zp*rv(5,i)/n_avg
    jy(x+1,y+1,z+1)=jy(x+1,y+1,z+1)+charge*xm*ym*zm*rv(5,i)/n_avg
    jz(x  ,y  ,z  )=jz(x  ,y  ,z  )+charge*xp*yp*zp*rv(6,i)/n_avg
    jz(x+1,y  ,z  )=jz(x+1,y  ,z  )+charge*xm*yp*zp*rv(6,i)/n_avg
    jz(x  ,y+1,z  )=jz(x  ,y+1,z  )+charge*xp*ym*zp*rv(6,i)/n_avg
    jz(x  ,y  ,z+1)=jz(x  ,y  ,z+1)+charge*xp*yp*zm*rv(6,i)/n_avg
    jz(x  ,y+1,z+1)=jz(x  ,y+1,z+1)+charge*xp*ym*zm*rv(6,i)/n_avg
    jz(x+1,y  ,z+1)=jz(x+1,y  ,z+1)+charge*xm*yp*zm*rv(6,i)/n_avg
    jz(x+1,y+1,z  )=jz(x+1,y+1,z  )+charge*xm*ym*zp*rv(6,i)/n_avg
    jz(x+1,y+1,z+1)=jz(x+1,y+1,z+1)+charge*xm*ym*zm*rv(6,i)/n_avg
    if ((mcx.ne.0) .or. (mcy.ne.0) .or. (mcz.ne.0)) then
      if ((my_pey==0 .and. y==0) .or. (my_pey==pey-1 .and. y==ny)) then
        tmp=yp; yp=ym; ym=tmp
        jx(x  ,y  ,z  )=jx(x  ,y  ,z  )+charge*xp*yp*zp*rv(4,i)/n_avg*mcx
        jx(x+1,y  ,z  )=jx(x+1,y  ,z  )+charge*xm*yp*zp*rv(4,i)/n_avg*mcx
        jx(x  ,y+1,z  )=jx(x  ,y+1,z  )+charge*xp*ym*zp*rv(4,i)/n_avg*mcx
        jx(x  ,y  ,z+1)=jx(x  ,y  ,z+1)+charge*xp*yp*zm*rv(4,i)/n_avg*mcx
        jx(x  ,y+1,z+1)=jx(x  ,y+1,z+1)+charge*xp*ym*zm*rv(4,i)/n_avg*mcx
        jx(x+1,y  ,z+1)=jx(x+1,y  ,z+1)+charge*xm*yp*zm*rv(4,i)/n_avg*mcx
        jx(x+1,y+1,z  )=jx(x+1,y+1,z  )+charge*xm*ym*zp*rv(4,i)/n_avg*mcx
        jx(x+1,y+1,z+1)=jx(x+1,y+1,z+1)+charge*xm*ym*zm*rv(4,i)/n_avg*mcx
        jy(x  ,y  ,z  )=jy(x  ,y  ,z  )+charge*xp*yp*zp*rv(5,i)/n_avg*mcy
        jy(x+1,y  ,z  )=jy(x+1,y  ,z  )+charge*xm*yp*zp*rv(5,i)/n_avg*mcy
        jy(x  ,y+1,z  )=jy(x  ,y+1,z  )+charge*xp*ym*zp*rv(5,i)/n_avg*mcy
        jy(x  ,y  ,z+1)=jy(x  ,y  ,z+1)+charge*xp*yp*zm*rv(5,i)/n_avg*mcy
        jy(x  ,y+1,z+1)=jy(x  ,y+1,z+1)+charge*xp*ym*zm*rv(5,i)/n_avg*mcy
        jy(x+1,y  ,z+1)=jy(x+1,y  ,z+1)+charge*xm*yp*zm*rv(5,i)/n_avg*mcy
        jy(x+1,y+1,z  )=jy(x+1,y+1,z  )+charge*xm*ym*zp*rv(5,i)/n_avg*mcy
        jy(x+1,y+1,z+1)=jy(x+1,y+1,z+1)+charge*xm*ym*zm*rv(5,i)/n_avg*mcy
        jz(x  ,y  ,z  )=jz(x  ,y  ,z  )+charge*xp*yp*zp*rv(6,i)/n_avg*mcz
        jz(x+1,y  ,z  )=jz(x+1,y  ,z  )+charge*xm*yp*zp*rv(6,i)/n_avg*mcz
        jz(x  ,y+1,z  )=jz(x  ,y+1,z  )+charge*xp*ym*zp*rv(6,i)/n_avg*mcz
        jz(x  ,y  ,z+1)=jz(x  ,y  ,z+1)+charge*xp*yp*zm*rv(6,i)/n_avg*mcz
        jz(x  ,y+1,z+1)=jz(x  ,y+1,z+1)+charge*xp*ym*zm*rv(6,i)/n_avg*mcz
        jz(x+1,y  ,z+1)=jz(x+1,y  ,z+1)+charge*xm*yp*zm*rv(6,i)/n_avg*mcz
        jz(x+1,y+1,z  )=jz(x+1,y+1,z  )+charge*xm*ym*zp*rv(6,i)/n_avg*mcz
        jz(x+1,y+1,z+1)=jz(x+1,y+1,z+1)+charge*xm*ym*zm*rv(6,i)/n_avg*mcz
      endif
    endif
  enddo
! transfer current in boundary layer to adjacent PEs
  call bound_rho_j(jx); call bound_rho_j(jy); call bound_rho_j(jz)
  call stop_timer(t_calc_j,n_calc_j)
End Subroutine calc_j

!-----------------------------------------------------------------------------
!                      Calculate Pressure Tensor
!-----------------------------------------------------------------------------

Subroutine calc_p(rv,np,m,charge,rho2,jx2,jy2,jz2,reset)
! calculates due to particle species with charge 'charge'
! rv vector with positions and velocities, np number of particles on PE
! current is added to jx,y,z(x,y,z), initialized to zero if reset=.true.

  Use pe_env
  Use partfield
  Use timer
  Implicit none
  Real(kind=prk), Dimension(6,maxparticles) :: rv
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1), intent(in) :: jx2,jy2,jz2,rho2
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez)
  Real(kind=drk) :: charge, xm, ym, zm, xp, yp, zp, tmp, m
  Integer :: np,i,x,y,z,mcx,mcy,mcz
  Logical :: reset
  
  call start_timer()
! reset arrays for current density
  if (reset) then
    do x=0,nx+1; do y=0,ny+1; do z=0,nz+1
      pxx(x,y,z)=0.; pyy(x,y,z)=0.; pzz(x,y,z)=0.
      pxz(x,y,z)=0.; pyz(x,y,z)=0.; pxy(x,y,z)=0.

    enddo; enddo; enddo
! pxx=0. ; pyy=0. ; pzz=0.
  else
! reset boundary layers only
    do x=0,nx+1; do y=0,ny+1
      pxx(x,y,0)=0.; pxx(x,y,nz+1)=0.
      pyy(x,y,0)=0.; pyy(x,y,nz+1)=0.
      pzz(x,y,0)=0.; pzz(x,y,nz+1)=0.
      pxz(x,y,0)=0.; pxz(x,y,nz+1)=0.
      pyz(x,y,0)=0.; pyz(x,y,nz+1)=0.
      pxy(x,y,0)=0.; pxy(x,y,nz+1)=0.
    enddo; enddo
    do x=0,nx+1; do z=0,nz+1
      pxx(x,0,z)=0.; pxx(x,ny+1,z)=0.
      pyy(x,0,z)=0.; pyy(x,ny+1,z)=0.
      pzz(x,0,z)=0.; pzz(x,ny+1,z)=0.
      pxz(x,0,z)=0.; pxz(x,ny+1,z)=0.
      pyz(x,0,z)=0.; pyz(x,ny+1,z)=0.
      pxy(x,0,z)=0.; pxy(x,ny+1,z)=0.
    enddo; enddo
    do y=0,ny+1; do z=0,nz+1
      pxx(0,y,z)=0.; pxx(nx+1,y,z)=0.
      pyy(0,y,z)=0.; pyy(nx+1,y,z)=0.
      pzz(0,y,z)=0.; pzz(nx+1,y,z)=0.
      pxz(0,y,z)=0.; pxz(nx+1,y,z)=0.
      pyz(0,y,z)=0.; pyz(nx+1,y,z)=0.
      pxy(0,y,z)=0.; pxy(nx+1,y,z)=0.
    enddo; enddo
  endif

! Virtual current is currently disabled. Code will only run with
! periodic boundary conditions.  
  call virtual_current(mcx,mcy,mcz)

! calculate pressure tensors
#if (nz == 1 && pez == 1) 
  z = 0; zm=1.; zp=0.
  do i=1, np
    x=int(rv(1,i)/dx-my_pex*nx+0.5)
    y=int(rv(2,i)/dy-my_pey*ny+0.5)
    xm=rv(1,i)/dx-my_pex*nx+0.5-x; xp=1.-xm
    ym=rv(2,i)/dy-my_pey*ny+0.5-y; yp=1.-ym
#else
  do i=1, np
    x=int(rv(1,i)/dx-my_pex*nx+0.5)
    y=int(rv(2,i)/dy-my_pey*ny+0.5)
    z=int(rv(3,i)/dz-my_pez*nz+0.5)
    xm=rv(1,i)/dx-my_pex*nx+0.5-x; xp=1.-xm
    ym=rv(2,i)/dy-my_pey*ny+0.5-y; yp=1.-ym
    zm=rv(3,i)/dz-my_pez*nz+0.5-z; zp=1.-zm
#endif
    pxx(x  ,y  ,z  )=pxx(x  ,y  ,z  )+m*xp*yp*zp*rv(4,i)*rv(4,i)/n_avg
    pxx(x+1,y  ,z  )=pxx(x+1,y  ,z  )+m*xm*yp*zp*rv(4,i)*rv(4,i)/n_avg
    pxx(x  ,y+1,z  )=pxx(x  ,y+1,z  )+m*xp*ym*zp*rv(4,i)*rv(4,i)/n_avg
    pxx(x  ,y  ,z+1)=pxx(x  ,y  ,z+1)+m*xp*yp*zm*rv(4,i)*rv(4,i)/n_avg
    pxx(x  ,y+1,z+1)=pxx(x  ,y+1,z+1)+m*xp*ym*zm*rv(4,i)*rv(4,i)/n_avg
    pxx(x+1,y  ,z+1)=pxx(x+1,y  ,z+1)+m*xm*yp*zm*rv(4,i)*rv(4,i)/n_avg
    pxx(x+1,y+1,z  )=pxx(x+1,y+1,z  )+m*xm*ym*zp*rv(4,i)*rv(4,i)/n_avg
    pxx(x+1,y+1,z+1)=pxx(x+1,y+1,z+1)+m*xm*ym*zm*rv(4,i)*rv(4,i)/n_avg
    pyy(x  ,y  ,z  )=pyy(x  ,y  ,z  )+m*xp*yp*zp*rv(5,i)*rv(5,i)/n_avg
    pyy(x+1,y  ,z  )=pyy(x+1,y  ,z  )+m*xm*yp*zp*rv(5,i)*rv(5,i)/n_avg
    pyy(x  ,y+1,z  )=pyy(x  ,y+1,z  )+m*xp*ym*zp*rv(5,i)*rv(5,i)/n_avg
    pyy(x  ,y  ,z+1)=pyy(x  ,y  ,z+1)+m*xp*yp*zm*rv(5,i)*rv(5,i)/n_avg
    pyy(x  ,y+1,z+1)=pyy(x  ,y+1,z+1)+m*xp*ym*zm*rv(5,i)*rv(5,i)/n_avg
    pyy(x+1,y  ,z+1)=pyy(x+1,y  ,z+1)+m*xm*yp*zm*rv(5,i)*rv(5,i)/n_avg
    pyy(x+1,y+1,z  )=pyy(x+1,y+1,z  )+m*xm*ym*zp*rv(5,i)*rv(5,i)/n_avg
    pyy(x+1,y+1,z+1)=pyy(x+1,y+1,z+1)+m*xm*ym*zm*rv(5,i)*rv(5,i)/n_avg
    pzz(x  ,y  ,z  )=pzz(x  ,y  ,z  )+m*xp*yp*zp*rv(6,i)*rv(6,i)/n_avg
    pzz(x+1,y  ,z  )=pzz(x+1,y  ,z  )+m*xm*yp*zp*rv(6,i)*rv(6,i)/n_avg
    pzz(x  ,y+1,z  )=pzz(x  ,y+1,z  )+m*xp*ym*zp*rv(6,i)*rv(6,i)/n_avg
    pzz(x  ,y  ,z+1)=pzz(x  ,y  ,z+1)+m*xp*yp*zm*rv(6,i)*rv(6,i)/n_avg
    pzz(x  ,y+1,z+1)=pzz(x  ,y+1,z+1)+m*xp*ym*zm*rv(6,i)*rv(6,i)/n_avg
    pzz(x+1,y  ,z+1)=pzz(x+1,y  ,z+1)+m*xm*yp*zm*rv(6,i)*rv(6,i)/n_avg
    pzz(x+1,y+1,z  )=pzz(x+1,y+1,z  )+m*xm*ym*zp*rv(6,i)*rv(6,i)/n_avg
    pzz(x+1,y+1,z+1)=pzz(x+1,y+1,z+1)+m*xm*ym*zm*rv(6,i)*rv(6,i)/n_avg
    pxz(x  ,y  ,z  )=pxz(x  ,y  ,z  )+m*xp*yp*zp*rv(4,i)*rv(6,i)/n_avg
    pxz(x+1,y  ,z  )=pxz(x+1,y  ,z  )+m*xm*yp*zp*rv(4,i)*rv(6,i)/n_avg
    pxz(x  ,y+1,z  )=pxz(x  ,y+1,z  )+m*xp*ym*zp*rv(4,i)*rv(6,i)/n_avg
    pxz(x  ,y  ,z+1)=pxz(x  ,y  ,z+1)+m*xp*yp*zm*rv(4,i)*rv(6,i)/n_avg
    pxz(x  ,y+1,z+1)=pxz(x  ,y+1,z+1)+m*xp*ym*zm*rv(4,i)*rv(6,i)/n_avg
    pxz(x+1,y  ,z+1)=pxz(x+1,y  ,z+1)+m*xm*yp*zm*rv(4,i)*rv(6,i)/n_avg
    pxz(x+1,y+1,z  )=pxz(x+1,y+1,z  )+m*xm*ym*zp*rv(4,i)*rv(6,i)/n_avg
    pxz(x+1,y+1,z+1)=pxz(x+1,y+1,z+1)+m*xm*ym*zm*rv(4,i)*rv(6,i)/n_avg
    pyz(x  ,y  ,z  )=pyz(x  ,y  ,z  )+m*xp*yp*zp*rv(5,i)*rv(6,i)/n_avg
    pyz(x+1,y  ,z  )=pyz(x+1,y  ,z  )+m*xm*yp*zp*rv(5,i)*rv(6,i)/n_avg
    pyz(x  ,y+1,z  )=pyz(x  ,y+1,z  )+m*xp*ym*zp*rv(5,i)*rv(6,i)/n_avg
    pyz(x  ,y  ,z+1)=pyz(x  ,y  ,z+1)+m*xp*yp*zm*rv(5,i)*rv(6,i)/n_avg
    pyz(x  ,y+1,z+1)=pyz(x  ,y+1,z+1)+m*xp*ym*zm*rv(5,i)*rv(6,i)/n_avg
    pyz(x+1,y  ,z+1)=pyz(x+1,y  ,z+1)+m*xm*yp*zm*rv(5,i)*rv(6,i)/n_avg
    pyz(x+1,y+1,z  )=pyz(x+1,y+1,z  )+m*xm*ym*zp*rv(5,i)*rv(6,i)/n_avg
    pyz(x+1,y+1,z+1)=pyz(x+1,y+1,z+1)+m*xm*ym*zm*rv(5,i)*rv(6,i)/n_avg
    pxy(x  ,y  ,z  )=pxy(x  ,y  ,z  )+m*xp*yp*zp*rv(4,i)*rv(5,i)/n_avg
    pxy(x+1,y  ,z  )=pxy(x+1,y  ,z  )+m*xm*yp*zp*rv(4,i)*rv(5,i)/n_avg
    pxy(x  ,y+1,z  )=pxy(x  ,y+1,z  )+m*xp*ym*zp*rv(4,i)*rv(5,i)/n_avg
    pxy(x  ,y  ,z+1)=pxy(x  ,y  ,z+1)+m*xp*yp*zm*rv(4,i)*rv(5,i)/n_avg
    pxy(x  ,y+1,z+1)=pxy(x  ,y+1,z+1)+m*xp*ym*zm*rv(4,i)*rv(5,i)/n_avg
    pxy(x+1,y  ,z+1)=pxy(x+1,y  ,z+1)+m*xm*yp*zm*rv(4,i)*rv(5,i)/n_avg
    pxy(x+1,y+1,z  )=pxy(x+1,y+1,z  )+m*xm*ym*zp*rv(4,i)*rv(5,i)/n_avg
    pxy(x+1,y+1,z+1)=pxy(x+1,y+1,z+1)+m*xm*ym*zm*rv(4,i)*rv(5,i)/n_avg
! Disable all but periodic boundary conditions.
    if ((mcx.ne.0) .or. (mcy.ne.0) .or. (mcz.ne.0)) then
      print *,'You should not get here because only periodic'
      print *,'boundary conditions are allowed'
      stop 789
      if ((my_pey==0 .and. y==0) .or. (my_pey==pey-1 .and. y==ny)) then
        tmp=yp; yp=ym; ym=tmp
        jx(x  ,y  ,z  )=jx(x  ,y  ,z  )+charge*xp*yp*zp*rv(4,i)/n_avg*mcx
        jx(x+1,y  ,z  )=jx(x+1,y  ,z  )+charge*xm*yp*zp*rv(4,i)/n_avg*mcx
        jx(x  ,y+1,z  )=jx(x  ,y+1,z  )+charge*xp*ym*zp*rv(4,i)/n_avg*mcx
        jx(x  ,y  ,z+1)=jx(x  ,y  ,z+1)+charge*xp*yp*zm*rv(4,i)/n_avg*mcx
        jx(x  ,y+1,z+1)=jx(x  ,y+1,z+1)+charge*xp*ym*zm*rv(4,i)/n_avg*mcx
        jx(x+1,y  ,z+1)=jx(x+1,y  ,z+1)+charge*xm*yp*zm*rv(4,i)/n_avg*mcx
        jx(x+1,y+1,z  )=jx(x+1,y+1,z  )+charge*xm*ym*zp*rv(4,i)/n_avg*mcx
        jx(x+1,y+1,z+1)=jx(x+1,y+1,z+1)+charge*xm*ym*zm*rv(4,i)/n_avg*mcx
        jy(x  ,y  ,z  )=jy(x  ,y  ,z  )+charge*xp*yp*zp*rv(5,i)/n_avg*mcy
        jy(x+1,y  ,z  )=jy(x+1,y  ,z  )+charge*xm*yp*zp*rv(5,i)/n_avg*mcy
        jy(x  ,y+1,z  )=jy(x  ,y+1,z  )+charge*xp*ym*zp*rv(5,i)/n_avg*mcy
        jy(x  ,y  ,z+1)=jy(x  ,y  ,z+1)+charge*xp*yp*zm*rv(5,i)/n_avg*mcy
        jy(x  ,y+1,z+1)=jy(x  ,y+1,z+1)+charge*xp*ym*zm*rv(5,i)/n_avg*mcy
        jy(x+1,y  ,z+1)=jy(x+1,y  ,z+1)+charge*xm*yp*zm*rv(5,i)/n_avg*mcy
        jy(x+1,y+1,z  )=jy(x+1,y+1,z  )+charge*xm*ym*zp*rv(5,i)/n_avg*mcy
        jy(x+1,y+1,z+1)=jy(x+1,y+1,z+1)+charge*xm*ym*zm*rv(5,i)/n_avg*mcy
        jz(x  ,y  ,z  )=jz(x  ,y  ,z  )+charge*xp*yp*zp*rv(6,i)/n_avg*mcz
        jz(x+1,y  ,z  )=jz(x+1,y  ,z  )+charge*xm*yp*zp*rv(6,i)/n_avg*mcz
        jz(x  ,y+1,z  )=jz(x  ,y+1,z  )+charge*xp*ym*zp*rv(6,i)/n_avg*mcz
        jz(x  ,y  ,z+1)=jz(x  ,y  ,z+1)+charge*xp*yp*zm*rv(6,i)/n_avg*mcz
        jz(x  ,y+1,z+1)=jz(x  ,y+1,z+1)+charge*xp*ym*zm*rv(6,i)/n_avg*mcz
        jz(x+1,y  ,z+1)=jz(x+1,y  ,z+1)+charge*xm*yp*zm*rv(6,i)/n_avg*mcz
        jz(x+1,y+1,z  )=jz(x+1,y+1,z  )+charge*xm*ym*zp*rv(6,i)/n_avg*mcz
        jz(x+1,y+1,z+1)=jz(x+1,y+1,z+1)+charge*xm*ym*zm*rv(6,i)/n_avg*mcz
      endif
    endif
  enddo
! transfer current in boundary layer to adjacent PEs
  call bound_rho_j(pxx); call bound_rho_j(pyy); call bound_rho_j(pzz)
  call bound_rho_j(pxz); call bound_rho_j(pyz); call bound_rho_j(pxy)

! Currently the pressures are actually pressure + convection term.
! Subtract out convection term.

  pxx = pxx - m*charge*jx2*jx2/rho2  !jx and rho have charge sign
  pyy = pyy - m*charge*jy2*jy2/rho2
  pzz = pzz - m*charge*jz2*jz2/rho2
  pxz = pxz - m*charge*jx2*jz2/rho2
  pyz = pyz - m*charge*jy2*jz2/rho2
  pxy = pxy - m*charge*jx2*jy2/rho2


  call stop_timer(t_calc_p,n_calc_p)
End Subroutine calc_p

!-----------------------------------------------------------------------------
!                         3D-boundary condition
!-----------------------------------------------------------------------------

! transfer charge or current density from boundary layer to adjacent box
! periodic boundary conditions, unless non_periodic() == .true. 
! (called by calc_rho and calc_j)

Subroutine bound_rho_j(f)
  Use pe_env ! Contains also mpif.h
  Implicit None

  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: f
  Real(kind=drk), Dimension(0:ny+1,0:nz+1) :: tmpxs, tmpxr
  Real(kind=drk), Dimension(1:nx,  0:nz+1) :: tmpys, tmpyr
  Real(kind=drk), Dimension(1:nx,  1:ny  ) :: tmpzs, tmpzr
  Logical :: non_periodic
  Integer :: x, y, z, icnt
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status
  Integer :: mpi_err


! delete exterior boundary layers in y
  if (non_periodic()) then
    if (my_pey==0) then
      do x=0,nx+1; do z=0,nz+1; f(x,0,z)=0.; enddo; enddo
    endif
    if (my_pey==pey-1) then
      do x=0,nx+1; do z=0,nz+1; f(x,ny+1,z)=0.; enddo; enddo
    endif
  endif

! Exchange in X-Direction

  icnt = (ny+2)*(nz+2)

  tmpxs(0:ny+1,0:nz+1) = f(0,0:ny+1,0:nz+1)
  call MPI_Sendrecv(tmpxs,icnt,mpi_drk,nb_pe(-1,0,0),111, &
                    tmpxr,icnt,mpi_drk,nb_pe( 1,0,0),111, &
                    MPI_COMM_WORLD,mpi_status,mpi_err)
  f(nx,0:ny+1,0:nz+1) = f(nx,0:ny+1,0:nz+1) + tmpxr(0:ny+1,0:nz+1)

  tmpxs(0:ny+1,0:nz+1) = f(nx+1,0:ny+1,0:nz+1)
  call MPI_Sendrecv(tmpxs,icnt,mpi_drk,nb_pe( 1,0,0),222, &
                    tmpxr,icnt,mpi_drk,nb_pe(-1,0,0),222, &
                    MPI_COMM_WORLD,mpi_status,mpi_err)
  f(1,0:ny+1,0:nz+1) = f(1,0:ny+1,0:nz+1) + tmpxr(0:ny+1,0:nz+1)

! Exchange in Y-Direction

  icnt = nx*(nz+2)

  tmpys(1:nx,0:nz+1) = f(1:nx,0,0:nz+1)
  call MPI_Sendrecv(tmpys,icnt,mpi_drk,nb_pe(0,-1,0),333, &
                    tmpyr,icnt,mpi_drk,nb_pe(0, 1,0),333, &
                    MPI_COMM_WORLD,mpi_status,mpi_err)
  f(1:nx,ny,0:nz+1) = f(1:nx,ny,0:nz+1) + tmpyr(1:nx,0:nz+1)

  tmpys(1:nx,0:nz+1) = f(1:nx,ny+1,0:nz+1)
  call MPI_Sendrecv(tmpys,icnt,mpi_drk,nb_pe(0, 1,0),444, &
                    tmpyr,icnt,mpi_drk,nb_pe(0,-1,0),444, &
                    MPI_COMM_WORLD,mpi_status,mpi_err)
  f(1:nx,1,0:nz+1) = f(1:nx,1,0:nz+1) + tmpyr(1:nx,0:nz+1)

! Exchange in Z-Direction

  icnt = nx*ny

  tmpzs(1:nx,1:ny) = f(1:nx,1:ny,0)
  call MPI_Sendrecv(tmpzs,icnt,mpi_drk,nb_pe(0,0,-1),555, &
                    tmpzr,icnt,mpi_drk,nb_pe(0,0, 1),555, &
                    MPI_COMM_WORLD,mpi_status,mpi_err)
  f(1:nx,1:ny,nz) = f(1:nx,1:ny,nz) + tmpzr(1:nx,1:ny)

  tmpzs(1:nx,1:ny) = f(1:nx,1:ny,nz+1)
  call MPI_Sendrecv(tmpzs,icnt,mpi_drk,nb_pe(0,0, 1),666, &
                    tmpzr,icnt,mpi_drk,nb_pe(0,0,-1),666, &
                    MPI_COMM_WORLD,mpi_status,mpi_err)
  f(1:nx,1:ny,1) = f(1:nx,1:ny,1) + tmpzr(1:nx,1:ny)

  return
End Subroutine bound_rho_j

!-----------------------------------------------------------------------------
!                      redistribute particles on PEs
!-----------------------------------------------------------------------------


Subroutine redistribute(rv,np)
! redistributes particles on PEs according to location in x space
! (called by stepx)

  Use partfield
  Use pe_env ! Contains also mpif.h
  Use timer
  Implicit none

  Real(kind=prk), Dimension(6,maxparticles) :: rv
  Integer(kind=prk), Dimension(maxparticles) :: rp

  Integer, Parameter :: maxbuff = 60000
  Real(kind=prk), Dimension(0:6,maxbuff) :: sendbuff, recvbuff
  Integer, Dimension(maxbuff) :: indarr

  Integer :: np, nps, npr, npnew, icnt, ip, i, iproc, iprocs, iprocr,j
  Integer, Dimension(0:pex*pey*pez-1) :: ipdist, ipsend, ipstart, ipsort
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status
  Integer :: mpi_err
  Integer :: maxi, indi, fillin
  Real(kind=drk) :: again, overflow

  Integer :: procval
  Logical :: slow_ex
  Real(kind=drk) :: bigjump

  Integer :: slow_index

!  call start_timer()

! NEW: Older versions of this procedure tried to pass particles
! between each possible pair of processors.  For large N that's painful.
! In many cases one need only consider nearest neighbor processors, the
! number of which is independent of the system size.

! There must be at least three processors in a dimension for the
! fast passing model to not pass the same particles twice.  That would
! be bad, so check for it here.  Otherwise set procval to be the number
! of nearest neighbors. Somewhat wasteful to check each run of the 
! subroutine, but easier.
  slow_ex = .false.
#if (nz ==1 && pez ==1)
  if (pex <3 .or. pey < 3) slow_ex = .true.
  procval = 8
#else
  if (pex <3 .or. pey < 3 .or. pez <3) slow_ex = .true.
  procval = 26
#endif

! Calculate new PE number for every particle
! MS: Commented out RJ's lines since if rp(i) is out of range, there's
! another problem that dumping the particles into a random processor
! will not solve.
  do i=1,np
#if (nz ==1 && pez ==1)
    rp(i)=int((rv(1,i)/lx)*pex)+&
	  int((rv(2,i)/ly)*pey)*pex
#else
    rp(i)=int((rv(1,i)/lx)*pex)+&
          int((rv(2,i)/ly)*pey)*pex+&
          int((rv(3,i)/lz)*pez)*pex*pey
#endif
!  RJ: Safety first !!!
!    if(rp(i) < 0) rp(i) = 0
!    if(rp(i) >= nprocs) rp(i) = nprocs-1
  enddo

  npr = 0  ! Number of particles in receive buffer, none at the beginning

  again = 1.

! MS: The start of the big loop.
  do

! Loop until redistribution is done
    if(again<=0.) exit

!   At this point, we have particles in the particle array rv
!   and (possibly) in the receive buffer recvbuff
!   Count how many particles go to which PE (in ipdist)
    ipdist(:) = 0

    do i=1,np
      ipdist(rp(i)) = ipdist(rp(i)) + 1
    enddo

    ! MS: If the relativistic flag is set it is impossible for particles
    ! to move more than one processor in a given timestep because their
    ! velocity is limited to c.  Make that explicit here, but make an
    ! exception for the case of t=0.

    ! MS: The short code only runs if relativistic and t > 0.  To make
    ! that work we define a dummy variable slow_index.  If it is greater
    ! than 0, run the slow version. 
#ifdef relativistic
    slow_index = 0
#else
    slow_index = 1
#endif
    if (t .eq. 0.) slow_index = slow_index + 1

    ! MS: Even if non-relativistic, the fast model is still possible.
    if (slow_index .gt. 0) then
       ! NEW.  If any particles are to be passed to a processor other than a
       ! neighbor, default to the old passing model.  This should be rare, 
       ! except at t = 0 (if the init routine has dumped particles on the 
       ! grid willy-nilly).
       bigjump = 0.
       do i = 0,nprocs-1
          if (.not. any(nb_pe==i)) then
             if (ipdist(i) .gt. 0) bigjump = bigjump+1.
             exit     ! No need to waste time once one is found.
          endif
       enddo
       if (totalsum(bigjump) > 0.) slow_ex = .true.
    endif


!   For every remote PE calculate how many particles for that PE go into
!   the send buffer and their start positions in send buffer.
!   If not all particles that must be sent fit into the send buffer,
!   the number of particles to be sent is limited to the size of the
!   send buffer in this turn and the again-flag is set.

    nps = 0        ! Number of particles in send buffer
    ipsend(:)  = 0 ! Number of particles to send to a given PE
    ipstart(:) = 0 ! Start position in send buffer (0 based!!!)
    again = 0.     ! Flag if we have to loop again

    do i=1,nprocs-1
      iproc = mod(myproc+i,nprocs)
      ipstart(iproc) = nps
      if(nps+ipdist(iproc) > maxbuff) then
!       Send buffer too small to exchange all particles in this turn
        ipsend(iproc) = maxbuff-nps
        again = 1.
      else
        ipsend(iproc) = ipdist(iproc)
      endif
      nps = nps+ipsend(iproc)
    enddo

!   Sort particles to be sent into the send buffer.
    ipsort(:) = 0  ! Auxiliary array for sorting.  
		   ! Equal to no. of particles going to a remote PE.
    npnew = 0      ! New number of particles in particle buffer

!   MS: Rewrote to improve efficiency.  Old version collapsed every
!   element of rv downwards as sorted items were removed.  indarr now
!   saves the indices of the particles moved to sendbuff.  The last
!   element of rv that was not itself sent is then moved to the spot
!   vacated by the last element of indarr.  And voila: A new rv that
!   contains particles from 1 to npnew.

!   New variable
    indi = 1
    do i=1,np
      ip = rp(i)
      ! MS: If particle belongs on this PE, or too many particles 
      ! have already been moved.
      if (ip == myproc .or. ipsort(ip)>=ipsend(ip)) then
        npnew = npnew+1
      else
        indarr(indi) = i
        indi = indi+1
        ipsort(ip) = ipsort(ip)+1
        sendbuff(0,  ipstart(ip)+ipsort(ip)) = ip
        sendbuff(1:6,ipstart(ip)+ipsort(ip)) = rv(1:6,i)
      endif
    enddo

!  WCM: Close gaps in rv left by particles which are in sendbuff now.
!       Vacant places are filled in ascending indarr-order by taking
!       particles from the end of rv.
    maxi = np
    indi=indi-1
    fillin=1
    if(indi > 0) then
      do
        if(maxi==indarr(indi)) then
	  maxi=maxi-1
	  indi=indi-1
	  if (indi<fillin) exit
        else
	  rv(1:6,indarr(fillin)) = rv(1:6,maxi)
	  rp(indarr(fillin)) =rp(maxi)
	  maxi=maxi-1
	  fillin = fillin+1
	  if (fillin > indi) exit
	endif
      enddo
    endif

    np = npnew

!   Exchange particles
! NEW.  If the slow exchange model is necessary, step over every processor.
    if (slow_ex) procval = nprocs-1

    do i=1,procval
! NEW.  Step over processors correctly.  For slow case, step pver
! everything but myself, for fast case step over nearest neighbors
! (defined in modules.F90)
      if (slow_ex) then
        iprocs = mod(myproc+i,nprocs)         ! Number of PE to send to
        iprocr = mod(myproc+nprocs-i,nprocs)  ! Number of PE to receive from
      else
#if (nz == 1 && pez == 1)
        iprocs = alt_nb_2d(i)
        iprocr = alt_nb_2d(9-i)
#else
        iprocs = alt_nb_3d(i)
        iprocr = alt_nb_3d(27-i)
#endif
      endif

      recvbuff(0,:) = -1.

! MS: Although it looks like only sendbuff(0,:) is being sent, that is
! actually just a pointer to the start of the array containing the
! 7*ipsend(iprocs) elements.  Also, 7*maxbuff is the maximum number of
! values that can be sent; in general there will be fewer.
      call MPI_Sendrecv ( &
        sendbuff(0,ipstart(iprocs)+1),7*ipsend(iprocs),mpi_prk,iprocs,i, &
        recvbuff,7*maxbuff,mpi_prk,iprocr,i,MPI_COMM_WORLD, &
        mpi_status,mpi_err)

! MS: Find the number of received particles. The first index of
! recvbuff is equal to the current processor number for number of
! particles received, -1 otherwise
      npr = 0
      do j=1,maxbuff
	if (recvbuff(0,j) .eq. -1.) exit
        npr = npr+1
      enddo

!   Check for particle buffer overflow, exit if that happens
      if (np+npr > maxparticles) then
        write(6,*) &
          '*****',myproc,' redistribute: particle buffer overflow, abort *****'
        call flush(6)
        call exitallpes()
      endif

!     Store the particles received in the particles buffer
      rv(1:6,np+1:np+npr) = recvbuff(1:6,1:npr)
      rp(np+1:np+npr)     = recvbuff(0,  1:npr)

      np = np+npr
      npr = 0

    enddo

    again = totalsum(again)
  enddo

! call stop_timer(t_redis,n_redis)

End Subroutine redistribute



!-----------------------------------------------------------------------------
!                      Sort particles for better cache usage
!-----------------------------------------------------------------------------

Subroutine partsort(rv,np)

  Use pe_env
  Use timer

  Implicit none
  Real(kind=prk), Dimension(6,maxparticles) :: rv
  Integer np

  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez)

  Integer(kind=4), Dimension(maxparticles) :: nc
  Real(kind=4), Dimension(6) :: rvx, rvs
  Integer :: x, y, z, i, lenold, lennew, istart, nl, ncx
  Integer, Dimension(0:(nx+1)*(ny+1)*(nz+1)) :: idx

  call start_timer()

! Calculate the cell number for each particle (in nc)

  do i=1,np
    x=rv(1,i)/dx-my_pex*nx+0.5
    y=rv(2,i)/dy-my_pey*ny+0.5
#if (nz ==1 && pez ==1) 
    nc(i) = x + y*(nx+1)
#else
    z=rv(3,i)/dz-my_pez*nz+0.5
    nc(i) = x + y*(nx+1) + z*(nx+1)*(ny+1)
#endif
  enddo

! Calculate how many particles go into every cell (in idx)

  idx(:) = 0
  do i=1,np
    idx(nc(i)) = idx(nc(i)) + 1
  enddo

! Calculate the starting index of every cell in the (sorted) particle array

  lenold = idx(0)
  idx(0) = 1
  do i=1,(nx+1)*(ny+1)*(nz+1)
    lennew = idx(i)
    idx(i) = idx(i-1)+lenold
    lenold = lennew
  enddo

! Sort the particles by cell number
! The tricky part here ist that we have to do that inplace.
! If we could use a new array for the sorted particles, the following would be
! as simple as:
!
!   do i=1,np
!     nl = idx(nc(i))
!     idx(nc(i)) = idx(nc(i))+1
!     rv_sorted(:,nl) = rv(:,i)
!   enddo
!
! The basic description of the inplace sorting is:
!
! Loop until all particles are sorted:
!   - fetch the next not yet sorted particle (and mark its location as empty)
!   - calculate the new location of that particle
!   - put the particle into its new location, fetch the particle which was
!     at that location and repeat the above unless new location was marked empty

  istart = 1
  do while(.true.)

    ! Search next particle which is not yet sorted

    do i=istart,np
      if(nc(i)>=0) exit
    enddo

    ! Check for end of sort, update istart

    if(i>np) exit
    istart = i

    ncx = nc(i)
    nc(i) = -1  ! Mark as empty
    rvx = rv(:,i)

    do while(ncx /= -1)

      ! Calculate the new (=sorted) location of current particle

      nl = idx(ncx)
      idx(ncx) = idx(ncx)+1

      ! Place current particle into its new location, make particle
      ! at this location to new current particle
      ! Loop will finish if the location was marked as empty

      ncx = nc(nl)
      nc(nl) = -2  ! Mark as sorted

      rvs(:) = rv(:,nl)
      rv(:,nl) = rvx(:)
      rvx(:) = rvs(:)

    enddo

  enddo

  call stop_timer(t_partsort,n_partsort)

End Subroutine partsort
