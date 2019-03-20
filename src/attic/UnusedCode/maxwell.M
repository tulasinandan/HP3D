!******************************************************************************
!                   3D particle code: Maxwell equations
!                           Andreas Zeiler, 1999
!
!                       LATEST CHANGE: May 3, 2001
!
!******************************************************************************

#include "param"

define(DELTA,`(&
  ($1(x+1,y,z)-2*$1(x,y,z)+$1(x-1,y,z))/dx**2+&
  ($1(x,y+1,z)-2*$1(x,y,z)+$1(x,y-1,z))/dy**2+&
  ($1(x,y,z+1)-2*$1(x,y,z)+$1(x,y,z-1))/dz**2&
  )')

define(DX2,`(&
  ($1(x+1,y,z)-$1(x-1,y,z))/(2*dx)&
  )')
define(DY2,`(&
  ($1(x,y+1,z)-$1(x,y-1,z))/(2*dy)&
  )')
define(DZ2,`(&
  ($1(x,y,z+1)-$1(x,y,z-1))/(2*dz)&
  )')

define(DX4,`(&
  (-$1(x+2,y,z)+8*$1(x+1,y,z)-8*$1(x-1,y,z)+$1(x-2,y,z))/(12*dx)&
  )')
define(DY4,`(&
  (-$1(x,y+2,z)+8*$1(x,y+1,z)-8*$1(x,y-1,z)+$1(x,y-2,z))/(12*dy)&
  )')
define(DZ4,`(&
  (-$1(x,y,z+2)+8*$1(x,y,z+1)-8*$1(x,y,z-1)+$1(x,y,z-2))/(12*dz)&
  )')

define(X,`
#ifdef fourth_order
    if (x .ne. 1 .and. x .ne. nx) then
      define(DX,defn(`DX4')) $1 undefine(`DX')
    else
      define(DX,defn(`DX2')) $1 undefine(`DX')
    endif
#else
    define(DX,defn(`DX2')) $1 undefine(`DX')
#endif
  ')
define(Y,`
#ifdef fourth_order
    if (y .ne. 1 .and. y .ne. ny) then
      define(DY,defn(`DY4')) $1 undefine(`DY')
    else
      define(DY,defn(`DY2')) $1 undefine(`DY')
    endif
#else
    define(DY,defn(`DY2')) $1 undefine(`DY')
#endif
  ')
define(Z,`
#ifdef fourth_order
    if (z .ne. 1 .and. z .ne. nz) then
      define(DZ,defn(`DZ4')) $1 undefine(`DZ')
    else
      define(DZ,defn(`DZ2')) $1 undefine(`DZ')
    endif
#else
    define(DZ,defn(`DZ2')) $1 undefine(`DZ')
#endif
  ')

!------------------------------------------------------------------------------
!                    step electric and magnetic field
!------------------------------------------------------------------------------

! advance E,B according to Maxwell's equations (leapfrog algorithm)
! e1, b1, j, rho: E, B at time t (on entry); at time t+dts (on exit)
! e2, b2, j, rho: E, B at time t-dts/2 (on entry)
!                 average of e1,b1 at time t and time t+dts (on exit)

Subroutine stepfield()
  Use pe_env
  Use partfield
  Use timer
  Implicit None

  Integer x,y,z
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez), &
                               dts=dt/(2*substeps), dtc=dts*c_2
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: bsx,bsy,bsz

  call start_timer()
! first half step: e2,b2 t-dts/2 --> t+dts/2
! Faraday's Law: dB/dt=-curl E
  call bound_e(e1x,e1y,e1z)
  do z=1,nz; do y=1,ny; do x=1,nx
    X(b2y(x,y,z)=b2y(x,y,z)+dts*DX(e1z); b2z(x,y,z)=b2z(x,y,z)-dts*DX(e1y))
    Y(b2x(x,y,z)=b2x(x,y,z)-dts*DY(e1z); b2z(x,y,z)=b2z(x,y,z)+dts*DY(e1x))
    Z(b2x(x,y,z)=b2x(x,y,z)+dts*DZ(e1y); b2y(x,y,z)=b2y(x,y,z)-dts*DZ(e1x))
  enddo; enddo; enddo
! Ampere's Law: dE/dt = curl B - j
  call bound_b(b1x,b1y,b1z)
  do z=1,nz; do y=1,ny; do x=1,nx
    e2x(x,y,z) = e2x(x,y,z) - dtc*jx(x,y,z)
    e2y(x,y,z) = e2y(x,y,z) - dtc*jy(x,y,z)
    e2z(x,y,z) = e2z(x,y,z) - dtc*jz(x,y,z)
    X(e2y(x,y,z)=e2y(x,y,z)-dtc*DX(b1z); e2z(x,y,z)=e2z(x,y,z)+dtc*DX(b1y))
    Y(e2x(x,y,z)=e2x(x,y,z)+dtc*DY(b1z); e2z(x,y,z)=e2z(x,y,z)-dtc*DY(b1x))
    Z(e2x(x,y,z)=e2x(x,y,z)-dtc*DZ(b1y); e2y(x,y,z)=e2y(x,y,z)+dtc*DZ(b1x))
  enddo; enddo; enddo

! second half step: e1,b1 t --> t+dts
! save b1x,b1y,b1z for connection of leapfrog grids
  do z=1,nz; do y=1,ny; do x=1,nx
    bsx(x,y,z) = b1x(x,y,z)
    bsy(x,y,z) = b1y(x,y,z)
    bsz(x,y,z) = b1z(x,y,z)
  enddo; enddo; enddo
! Faraday's Law: dB/dt=-curl E
  call bound_e(e2x,e2y,e2z)
  do z=1,nz; do y=1,ny; do x=1,nx
    X(b1y(x,y,z)=b1y(x,y,z)+dts*DX(e2z); b1z(x,y,z)=b1z(x,y,z)-dts*DX(e2y))
    Y(b1x(x,y,z)=b1x(x,y,z)-dts*DY(e2z); b1z(x,y,z)=b1z(x,y,z)+dts*DY(e2x))
    Z(b1x(x,y,z)=b1x(x,y,z)+dts*DZ(e2y); b1y(x,y,z)=b1y(x,y,z)-dts*DZ(e2x))
  enddo; enddo; enddo
! save e1x,e1y,e1z for connection of leapfrog grids
  do z=1,nz; do y=1,ny; do x=1,nx
    e2x(x,y,z) = e1x(x,y,z)
    e2y(x,y,z) = e1y(x,y,z)
    e2z(x,y,z) = e1z(x,y,z)
  enddo; enddo; enddo
! Ampere's Law: dE/dt = curl B - j
  call bound_b(b2x,b2y,b2z)
  do z=1,nz; do y=1,ny; do x=1,nx
    e1x(x,y,z) = e1x(x,y,z) - dtc*jx(x,y,z)
    e1y(x,y,z) = e1y(x,y,z) - dtc*jy(x,y,z)
    e1z(x,y,z) = e1z(x,y,z) - dtc*jz(x,y,z)
    X(e1y(x,y,z)=e1y(x,y,z)-dtc*DX(b2z); e1z(x,y,z)=e1z(x,y,z)+dtc*DX(b2y))
    Y(e1x(x,y,z)=e1x(x,y,z)+dtc*DY(b2z); e1z(x,y,z)=e1z(x,y,z)-dtc*DY(b2x))
    Z(e1x(x,y,z)=e1x(x,y,z)-dtc*DZ(b2y); e1y(x,y,z)=e1y(x,y,z)+dtc*DZ(b2x))
  enddo; enddo; enddo
! connect leapfrog grids
  do z=1,nz; do y=1,ny; do x=1,nx
    e2x(x,y,z) = (e2x(x,y,z)+e1x(x,y,z))/2
    e2y(x,y,z) = (e2y(x,y,z)+e1y(x,y,z))/2
    e2z(x,y,z) = (e2z(x,y,z)+e1z(x,y,z))/2
    b2x(x,y,z) = (bsx(x,y,z)+b1x(x,y,z))/2
    b2y(x,y,z) = (bsy(x,y,z)+b1y(x,y,z))/2
    b2z(x,y,z) = (bsz(x,y,z)+b1z(x,y,z))/2
  enddo; enddo; enddo

  call stop_timer(t_stepfield,n_stepfield)
End Subroutine stepfield

!------------------------------------------------------------------------------
!                      fix divergence of E and B
!------------------------------------------------------------------------------

! correct E to satisfy Poisson's equation and B to be divergence free
! B is treated only if div_b=.true.

Subroutine divergence(ex,ey,ez,bx,by,bz,div_b)
  Use pe_env
  Use partfield
  Use timer
  Implicit None

  Integer x,y,z
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez)
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: ex,ey,ez,bx,by,bz
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: phi, rho_c
  Logical :: div_b

  call start_timer()

  phi = 0.

! if div_b make B divergence free (same algorithm as for Poisson equation)
  if (div_b) then
    call bound_b(bx,by,bz)
    do z=1,nz; do y=1,ny; do x=1,nx
      rho_c(x,y,z) = -DX2(bx)-DY2(by)-DZ2(bz)
    enddo; enddo; enddo
    call poisson(phi,rho_c,real(eps,drk))
    call bound_phi(0,0,phi)
    do z=1,nz; do y=1,ny; do x=1,nx
      X(bx(x,y,z) = bx(x,y,z) + DX(phi))
      Y(by(x,y,z) = by(x,y,z) + DY(phi))
      Z(bz(x,y,z) = bz(x,y,z) + DZ(phi))
    enddo; enddo; enddo
  endif

! Poisson's equation
! j and rho do not exactly satisfy the continuity equation:
! solve Poisson's equation   div grad phi = rho_c, rho_c = rho - div E
! correct E by adding grad phi
  call bound_e(ex,ey,ez)
  do z=1,nz; do y=1,ny; do x=1,nx
    rho_c(x,y,z) = c_2*rho(x,y,z) - (DX2(ex)+DY2(ey)+DZ2(ez))
  enddo; enddo; enddo
  call poisson(phi,rho_c,real(eps,drk))
! update E
  call bound_phi(0,0,phi)
  do z=1,nz; do y=1,ny; do x=1,nx
    X(ex(x,y,z) = ex(x,y,z) + DX(phi))
    Y(ey(x,y,z) = ey(x,y,z) + DY(phi))
    Z(ez(x,y,z) = ez(x,y,z) + DZ(phi))
  enddo; enddo; enddo
  call stop_timer(t_divergence,n_divergence)
End Subroutine divergence



