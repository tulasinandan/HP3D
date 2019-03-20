!******************************************************************************
!          3D particle code: fluid equations for hybrid version
!                           Andreas Zeiler, 2001
!
!                       LATEST CHANGE: July 20, 2001
!
!******************************************************************************

#include "param"

#ifndef D_b
#define D_b 0.
#endif
#ifndef D_pe
#define D_pe 0.
#endif
#ifndef n_base
#define n_base 0.
#endif

define(DELTA,`(&
  ($1(x+1,y,z)-2*$1(x,y,z)+$1(x-1,y,z))/dx**2+&
  ($1(x,y+1,z)-2*$1(x,y,z)+$1(x,y-1,z))/dy**2+&
  ($1(x,y,z+1)-2*$1(x,y,z)+$1(x,y,z-1))/dz**2&
  )')

define(HELMHOLTZ,`(&
  $1(x,y,z)-d_e2*(&
  ($1(x+1,y,z)-2*$1(x,y,z)+$1(x-1,y,z))/dx**2+&
  ($1(x,y+1,z)-2*$1(x,y,z)+$1(x,y-1,z))/dy**2+&
  ($1(x,y,z+1)-2*$1(x,y,z)+$1(x,y,z-1))/dz**2&
  ) )')

define(DX,`(&
  ($1(x+1,y,z)-$1(x-1,y,z))/(2*dx)&
  )')
define(DY,`(&
  ($1(x,y+1,z)-$1(x,y-1,z))/(2*dy)&
  )')
define(DZ,`(&
  ($1(x,y,z+1)-$1(x,y,z-1))/(2*dz)&
  )')

define(VECX,`(&
  $1y(x,y,z)*$2z(x,y,z)-$1z(x,y,z)*$2y(x,y,z)&
  )')
define(VECY,`(&
  $1z(x,y,z)*$2x(x,y,z)-$1x(x,y,z)*$2z(x,y,z)&
  )')
define(VECZ,`(&
  $1x(x,y,z)*$2y(x,y,z)-$1y(x,y,z)*$2x(x,y,z)&
  )')

!------------------------------------------------------------------------------
!                step magnetic field and electron pressure
!------------------------------------------------------------------------------

! advance B and pe according to fluid equations (leapfrog algorithm)
! b1, pe1: B, pe at time t (on entry); at time t+dts (on exit)
! b2, pe2: B, pe at time t-dts/2 (on entry)
!                 average of b1, pe1 at time t and time t+dts (on exit)

Subroutine stephybrid(fixed_cycles)
  Use pe_env
  Use partfield
  Use timer
  Implicit None

  Integer x,y,z,i,dampcycles
  Logical :: fixed_cycles
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), & 
                               dz=lz/(nz*pez), dts=dt/(2*substeps)
!     jt: temporary array to store the total current j = curl B
  Real(kind=drk), Save, Dimension(0:nx+1,0:ny+1,0:nz+1) :: jtx,jty,jtz

  call start_timer()
#ifdef isothermal
  do z=1,nz; do y=1,ny; do x=1,nx
    pe1(x,y,z)=rho(x,y,z)*T_e
    pe2(x,y,z)=rho(x,y,z)*T_e
  enddo; enddo; enddo
#endif
! first half step: b2,pe2 t-dts/2 --> t+dts/2
!     calculate E = ((curl b)x(1-d_e^2 delt^2)B - (ji x B) - delt pe)/n
!     esx,esy,esz are used as temporary arrays to store B'
  call bound_b(b1x,b1y,b1z)
  call bound_pe(pe1)
  do z=1,nz; do y=1,ny; do x=1,nx
    esx(x,y,z)=HELMHOLTZ(b1x) 
    esy(x,y,z)=HELMHOLTZ(b1y) 
    esz(x,y,z)=HELMHOLTZ(b1z)
    jtx(x,y,z)=DY(b1z)-DZ(b1y) 
    jty(x,y,z)=DZ(b1x)-DX(b1z) 
    jtz(x,y,z)=DX(b1y)-DY(b1x)
    e1x(x,y,z)=(VECX(jt,es)-VECX(j,b1)-DX(pe1))/rho(x,y,z)
    e1y(x,y,z)=(VECY(jt,es)-VECY(j,b1)-DY(pe1))/rho(x,y,z)
    e1z(x,y,z)=(VECZ(jt,es)-VECZ(j,b1)-DZ(pe1))/rho(x,y,z)
  enddo; enddo; enddo
!     calculate B2' from B2, advance B2'
  call bound_e(e1x,e1y,e1z)
  call bound_b(b2x,b2y,b2z)
  do z=1,nz; do y=1,ny; do x=1,nx
    esx(x,y,z)=HELMHOLTZ(b2x); esx(x,y,z)=esx(x,y,z)-dts*(DY(e1z)-DZ(e1y))
    esy(x,y,z)=HELMHOLTZ(b2y); esy(x,y,z)=esy(x,y,z)-dts*(DZ(e1x)-DX(e1z))
    esz(x,y,z)=HELMHOLTZ(b2z); esz(x,y,z)=esz(x,y,z)-dts*(DX(e1y)-DY(e1x))
  enddo; enddo; enddo
!     invert Helmholtz equation to obtain B2
  call helmholtz(b2x,esx,real(eps,drk),fixed_cycles) 
  call helmholtz(b2y,esy,real(eps,drk),fixed_cycles)
  call helmholtz(b2z,esz,real(eps,drk),fixed_cycles)
#ifndef isothermal
!     calculate electron velocity, multiply by rho**gam
  do z=1,nz; do y=1,ny; do x=1,nx
    jtx(x,y,z)=(jx(x,y,z)-jtx(x,y,z))*rho(x,y,z)**(gam-1)
    jty(x,y,z)=(jy(x,y,z)-jty(x,y,z))*rho(x,y,z)**(gam-1)
    jtz(x,y,z)=(jz(x,y,z)-jtz(x,y,z))*rho(x,y,z)**(gam-1)
  enddo; enddo; enddo
!     advance pe2
  do z=1,nz; do y=1,ny; do x=1,nx
    pe1(x,y,z)=pe1(x,y,z)*rho(x,y,z)**(-gam)
  enddo; enddo; enddo
  call bound_pe(pe1)
  do z=1,nz; do y=1,ny; do x=1,nx
    pe2(x,y,z)=pe2(x,y,z)+dts*( &
        jtx(x,y,z)*DX(pe1)+jty(x,y,z)*DY(pe1)+jtz(x,y,z)*DZ(pe1))
  enddo; enddo; enddo
  do z=1,nz; do y=1,ny; do x=1,nx
    pe1(x,y,z)=pe1(x,y,z)*rho(x,y,z)**(gam)
  enddo; enddo; enddo
#endif
! fourth order hyperviscosity for B and p_e
  dampcycles=int(max(D_b,D_pe)*dts*(3.142/min(dx,dy,dz))**4*20.+1.)
  do i=1,dampcycles
! use jtx, jty, jtz as temporary arrays
    call bound_b(b2x,b2y,b2z)
    do z=1,nz; do y=1,ny; do x=1,nx
      jtx(x,y,z)=DELTA(b2x) 
      jty(x,y,z)=DELTA(b2y) 
      jtz(x,y,z)=DELTA(b2z)
    enddo; enddo; enddo
    call bound_b(jtx,jty,jtz)
    do z=1,nz; do y=1,ny; do x=1,nx
      b2x(x,y,z)=b2x(x,y,z)-D_b*dts/dampcycles*DELTA(jtx) 
      b2y(x,y,z)=b2y(x,y,z)-D_b*dts/dampcycles*DELTA(jty) 
      b2z(x,y,z)=b2z(x,y,z)-D_b*dts/dampcycles*DELTA(jtz)
    enddo; enddo; enddo
#ifndef isothermal
    call bound_pe(pe2)
    do z=1,nz; do y=1,ny; do x=1,nx
      jtx(x,y,z)=DELTA(pe2) 
    enddo; enddo; enddo
    call bound_pe(jtx)
    do z=1,nz; do y=1,ny; do x=1,nx
      pe2(x,y,z)=pe2(x,y,z)-D_pe*dts/dampcycles*DELTA(jtx) 
    enddo; enddo; enddo
#endif
  enddo
! second half step: b1,pe1 t --> t+dts
!     calculate E = ((curl b)x(1-d_e^2 delt^2)B - (ji x B) - delt pe)/n
  call bound_b(b2x,b2y,b2z)
  call bound_pe(pe2)
  do z=1,nz; do y=1,ny; do x=1,nx
    esx(x,y,z)=HELMHOLTZ(b2x) 
    esy(x,y,z)=HELMHOLTZ(b2y) 
    esz(x,y,z)=HELMHOLTZ(b2z)
    jtx(x,y,z)=DY(b2z)-DZ(b2y) 
    jty(x,y,z)=DZ(b2x)-DX(b2z) 
    jtz(x,y,z)=DX(b2y)-DY(b2x)
    e1x(x,y,z)=(VECX(jt,es)-VECX(j,b2)-DX(pe2))/rho(x,y,z)
    e1y(x,y,z)=(VECY(jt,es)-VECY(j,b2)-DY(pe2))/rho(x,y,z)
    e1z(x,y,z)=(VECZ(jt,es)-VECZ(j,b2)-DZ(pe2))/rho(x,y,z)
  enddo; enddo; enddo
!     calculate B1' from B1, advance B1'
  call bound_e(e1x,e1y,e1z)
  do z=1,nz; do y=1,ny; do x=1,nx
    esx(x,y,z)=HELMHOLTZ(b1x); esx(x,y,z)=esx(x,y,z)-dts*(DY(e1z)-DZ(e1y))
    esy(x,y,z)=HELMHOLTZ(b1y); esy(x,y,z)=esy(x,y,z)-dts*(DZ(e1x)-DX(e1z))
    esz(x,y,z)=HELMHOLTZ(b1z); esz(x,y,z)=esz(x,y,z)-dts*(DX(e1y)-DY(e1x))
  enddo; enddo; enddo
!     save B1 to connect leapfrog grids
  do z=1,nz; do y=1,ny; do x=1,nx
    b2x(x,y,z)=b1x(x,y,z); b2y(x,y,z)=b1y(x,y,z); b2z(x,y,z)=b1z(x,y,z)
  enddo; enddo; enddo
!     invert Helmholtz equation to obtain B1
  call helmholtz(b1x,esx,real(eps,drk),fixed_cycles) 
  call helmholtz(b1y,esy,real(eps,drk),fixed_cycles)
  call helmholtz(b1z,esz,real(eps,drk),fixed_cycles)
#ifndef isothermal
!     calculate electron velocity, multiply by rho**gam
  do z=1,nz; do y=1,ny; do x=1,nx
    jtx(x,y,z)=(jx(x,y,z)-jtx(x,y,z))*rho(x,y,z)**(gam-1)
    jty(x,y,z)=(jy(x,y,z)-jty(x,y,z))*rho(x,y,z)**(gam-1)
    jtz(x,y,z)=(jz(x,y,z)-jtz(x,y,z))*rho(x,y,z)**(gam-1)
  enddo; enddo; enddo
!     save pe1 to connect grids, advance pe1 
  do z=1,nz; do y=1,ny; do x=1,nx
    pe2(x,y,z)=pe2(x,y,z)*rho(x,y,z)**(-gam)
  enddo; enddo; enddo
  call bound_pe(pe2)
  do z=1,nz; do y=1,ny; do x=1,nx
    esx(x,y,z)=pe1(x,y,z)
    pe1(x,y,z)=pe1(x,y,z)+dts*( &
        jtx(x,y,z)*DX(pe2)+jty(x,y,z)*DY(pe2)+jtz(x,y,z)*DZ(pe2))
  enddo; enddo; enddo
  do z=1,nz; do y=1,ny; do x=1,nx
    pe2(x,y,z)=pe2(x,y,z)*rho(x,y,z)**(gam)
  enddo; enddo; enddo
#endif
! fourth order hyperviscosity for B and p_e
  dampcycles=int(max(D_b,D_pe)*dts*(3.142/min(dx,dy,dz))**4*20.+1.)
  do i=1,dampcycles
    call bound_b(b1x,b1y,b1z)
    do z=1,nz; do y=1,ny; do x=1,nx
      jtx(x,y,z)=DELTA(b1x)
      jty(x,y,z)=DELTA(b1y)
      jtz(x,y,z)=DELTA(b1z)
    enddo; enddo; enddo
    call bound_b(jtx,jty,jtz)
    do z=1,nz; do y=1,ny; do x=1,nx
      b1x(x,y,z)=b1x(x,y,z)-D_b*dts/dampcycles*DELTA(jtx)
      b1y(x,y,z)=b1y(x,y,z)-D_b*dts/dampcycles*DELTA(jty)
      b1z(x,y,z)=b1z(x,y,z)-D_b*dts/dampcycles*DELTA(jtz)
    enddo; enddo; enddo
#ifndef isothermal
    call bound_pe(pe1)
    do z=1,nz; do y=1,ny; do x=1,nx
      jtx(x,y,z)=DELTA(pe1)
    enddo; enddo; enddo
    call bound_pe(jtx)
    do z=1,nz; do y=1,ny; do x=1,nx
      pe1(x,y,z)=pe1(x,y,z)-D_pe*dts/dampcycles*DELTA(jtx)
    enddo; enddo; enddo
#endif
  enddo
! connect leapfrog grids
  do z=1,nz; do y=1,ny; do x=1,nx
    b2x(x,y,z) = (b2x(x,y,z)+b1x(x,y,z))/2
    b2y(x,y,z) = (b2y(x,y,z)+b1y(x,y,z))/2
    b2z(x,y,z) = (b2z(x,y,z)+b1z(x,y,z))/2
#ifndef isothermal
    pe2(x,y,z) = (esx(x,y,z)+pe1(x,y,z))/2
#endif
  enddo; enddo; enddo

  call stop_timer(t_stephybrid,n_stephybrid)
End Subroutine stephybrid

!------------------------------------------------------------------------------
!                      fix divergence of B
!------------------------------------------------------------------------------

! correct B to be divergence free

Subroutine divergence(bx,by,bz)
  Use pe_env
  Use partfield
  Implicit None

  Integer x,y,z
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez)
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: bx,by,bz

! use esx and esy as temporary arrays
  esx(:,:,:) = 0

  call bound_b(bx,by,bz)
  do z=1,nz; do y=1,ny; do x=1,nx
    esy(x,y,z) = - (DX(bx)+DY(by)+DZ(bz))
  enddo; enddo; enddo
  call poisson(esx,esy,real(eps,drk))
  call bound_phi(0,0,esx)
  do z=1,nz; do y=1,ny; do x=1,nx
    bx(x,y,z) = bx(x,y,z) + DX(esx)
    by(x,y,z) = by(x,y,z) + DY(esx)
    bz(x,y,z) = bz(x,y,z) + DZ(esx)
  enddo; enddo; enddo
End Subroutine divergence

!------------------------------------------------------------------------------
!                        particle electric field
!------------------------------------------------------------------------------

! calculate electric field for particle stepper
! E'' = (j_e x B - delt p_e)/n, result in e1x,e1y,e1z

Subroutine calc_e()
  Use pe_env
  Use partfield
  Use timer
  Implicit None

  Integer x,y,z
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), & 
                               dz=lz/(nz*pez), dts=dt/(2*substeps)

  call start_timer()
  call bound_b(b1x,b1y,b1z)
  call bound_pe(pe1)
! use esx,esy,esz as temporary arrays to store the electron current
  do z=1,nz; do y=1,ny; do x=1,nx
    esx(x,y,z)=DY(b1z)-DZ(b1y)-jx(x,y,z)
    esy(x,y,z)=DZ(b1x)-DX(b1z)-jy(x,y,z)
    esz(x,y,z)=DX(b1y)-DY(b1x)-jz(x,y,z)
    e1x(x,y,z)=(VECX(es,b1)-DX(pe1))/rho(x,y,z)
    e1y(x,y,z)=(VECY(es,b1)-DY(pe1))/rho(x,y,z)
    e1z(x,y,z)=(VECZ(es,b1)-DZ(pe1))/rho(x,y,z)
  enddo; enddo; enddo
  call stop_timer(t_calc_e,n_calc_e)
End Subroutine calc_e

!------------------------------------------------------------------------------
!                               base density
!------------------------------------------------------------------------------

Subroutine add_base_density()
! add uniform density n_base to avoid division by zero in hybrid code

  Use pe_env
  Use partfield
  Implicit none
  Integer :: x,y,z

  do z=1,nz; do y=1,ny; do x=1,nx
    rho(x,y,z)=rho(x,y,z)+n_base
  enddo; enddo; enddo
End Subroutine add_base_density
