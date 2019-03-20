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
#ifndef mu
#define mu 0.
#endif
#ifndef n_base
#define n_base 0.
#endif

#ifndef EXPANSION
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

define(HX,`(&
  vhx(x,y,z)*($1(x+1,y,z)-$1(x,y,z))&
  -vhx(x-1,y,z)*($1(x,y,z)-$1(x-1,y,z))&
  )')

define(HY,`(&
  vhy(x,y,z)*($1(x,y+1,z)-$1(x,y,z))&
  -vhy(x,y-1,z)*($1(x,y,z)-$1(x,y-1,z))&
  )')

define(HZ,`(&
  vhz(x,y,z)*($1(x,y,z+1)-$1(x,y,z))&
  -vhz(x,y,z-1)*($1(x,y,z)-$1(x,y,z-1))&
  )')
#else
define(DELTA,`(&
  ($1(x+1,y,z)-2*$1(x,y,z)+$1(x-1,y,z))/(ww**2*dx**2)+&
  ($1(x,y+1,z)-2*$1(x,y,z)+$1(x,y-1,z))/(aa**2*dy**2)+&
  ($1(x,y,z+1)-2*$1(x,y,z)+$1(x,y,z-1))/(aa**2*dz**2)&
  )')

define(HELMHOLTZ,`(&
  $1(x,y,z)-d_e2*(&
  ($1(x+1,y,z)-2*$1(x,y,z)+$1(x-1,y,z))/(ww**2*dx**2)+&
  ($1(x,y+1,z)-2*$1(x,y,z)+$1(x,y-1,z))/(aa**2*dy**2)+&
  ($1(x,y,z+1)-2*$1(x,y,z)+$1(x,y,z-1))/(aa**2*dz**2)&
  ) )')

define(DX,`(&
  ($1(x+1,y,z)-$1(x-1,y,z))/(ww*2*dx)&
  )')
define(DY,`(&
  ($1(x,y+1,z)-$1(x,y-1,z))/(aa*2*dy)&
  )')
define(DZ,`(&
  ($1(x,y,z+1)-$1(x,y,z-1))/(aa*2*dz)&
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

define(HX,`(&
  vhx(x,y,z)*($1(x+1,y,z)-$1(x,y,z))&
  -vhx(x-1,y,z)*($1(x,y,z)-$1(x-1,y,z))&
  )')

define(HY,`(&
  vhy(x,y,z)*($1(x,y+1,z)-$1(x,y,z))&
  -vhy(x,y-1,z)*($1(x,y,z)-$1(x,y-1,z))&
  )')

define(HZ,`(&
  vhz(x,y,z)*($1(x,y,z+1)-$1(x,y,z))&
  -vhz(x,y,z-1)*($1(x,y,z)-$1(x,y,z-1))&
  )')
#endif

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
#ifdef EXPANSION
  Use expansion_params
#endif
  Implicit None

  Integer x,y,z,i,dampcycles
  Logical :: fixed_cycles
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), & 
                               dz=lz/(nz*pez), dts=dt/(2*substeps)
!     jt: temporary array to store the total current j = curl B
  Real(kind=drk), Save, Dimension(0:nx+1,0:ny+1,0:nz+1) :: jtx,jty,jtz,&
                                                           jt1x,jt1y,jt1z,&
                                                           jt2x,jt2y,jt2z
!     vh: temporary arrays to store (constant * average velocity) part 
!         of velocity-dependent hyperviscosity
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: vhx,vhy,vhz,vex,vey,vez,te

!---- aleph here is a multiplicative factor to turn off the electron pressure 
!---- in the electric field if need be.
#ifdef NO_ELEC_PHYS
  Real(kind=drk), Parameter :: aleph=0._8
#else
  Real(kind=drk), Parameter :: aleph=1._8 
#endif

  call start_timer()
#ifdef NO_ELEC_PHYS
#else
#ifdef isothermal
  do z=1,nz; do y=1,ny; do x=1,nx
    pe1(x,y,z)=rho(x,y,z)*T_e
    pe2(x,y,z)=rho(x,y,z)*T_e
  enddo; enddo; enddo
#endif
#ifdef barotropic
  do z=1,nz; do y=1,ny; do x=1,nx
    pe1(x,y,z)=T_e*rho(x,y,z)**gam
    pe2(x,y,z)=T_e*rho(x,y,z)**gam
  enddo; enddo; enddo
#endif
  call bound_pe(pe1)
#endif

  call bound_b(b1x,b1y,b1z)
!Calculate J_tot
  do z=1,nz; do y=1,ny; do x=1,nx
    jtx(x,y,z)=DY(b1z)-DZ(b1y) 
    jty(x,y,z)=DZ(b1x)-DX(b1z) 
    jtz(x,y,z)=DX(b1y)-DY(b1x)
  enddo; enddo; enddo

! first half step: b2,pe2 t-dts/2 --> t+dts/2
!     calculate E = ((curl b)x(1-d_e^2 delt^2)B - (ji x B) - delt pe)/n
!     esx,esy,esz are used as temporary arrays to store Bprime
  do z=1,nz; do y=1,ny; do x=1,nx
    esx(x,y,z)=HELMHOLTZ(b1x) 
    esy(x,y,z)=HELMHOLTZ(b1y) 
    esz(x,y,z)=HELMHOLTZ(b1z)
    e1x(x,y,z)=(VECX(jt,es)-VECX(j,b1))/rho(x,y,z) 
    e1y(x,y,z)=(VECY(jt,es)-VECY(j,b1))/rho(x,y,z) 
    e1z(x,y,z)=(VECZ(jt,es)-VECZ(j,b1))/rho(x,y,z) 
!   e1x(x,y,z)=(V_ECX(jt,es)-V_ECX(j,b1)-aleph*D_X(pe1))/rho(x,y,z) 
!   e1y(x,y,z)=(V_ECY(jt,es)-V_ECY(j,b1)-aleph*D_Y(pe1))/rho(x,y,z) 
!   e1z(x,y,z)=(V_ECZ(jt,es)-V_ECZ(j,b1)-aleph*D_Z(pe1))/rho(x,y,z) 
  enddo; enddo; enddo

!     calculate B2' from B2, advance B2'
  call bound_e(e1x,e1y,e1z)
  call bound_b(b2x,b2y,b2z)
#ifdef DRIVE_B
    if (myproc == 0) write(6,*) 'cos(4wt)      ',cos(FOUROMEGA*t)
#endif
  do z=1,nz; do y=1,ny; do x=1,nx
#ifdef DRIVE_B
    esx(x,y,z)=HELMHOLTZ(b2x); esx(x,y,z)=esx(x,y,z)-dts*(DY(e1z)-DZ(e1y)-cos(FOUROMEGA*t)*dr_bx(x,y,z))
    esy(x,y,z)=HELMHOLTZ(b2y); esy(x,y,z)=esy(x,y,z)-dts*(DZ(e1x)-DX(e1z)-cos(FOUROMEGA*t)*dr_by(x,y,z))
    esz(x,y,z)=HELMHOLTZ(b2z); esz(x,y,z)=esz(x,y,z)-dts*(DX(e1y)-DY(e1x)-cos(FOUROMEGA*t)*dr_bz(x,y,z))

#else
    esx(x,y,z)=HELMHOLTZ(b2x); esx(x,y,z)=esx(x,y,z)-dts*(DY(e1z)-DZ(e1y))
    esy(x,y,z)=HELMHOLTZ(b2y); esy(x,y,z)=esy(x,y,z)-dts*(DZ(e1x)-DX(e1z))
    esz(x,y,z)=HELMHOLTZ(b2z); esz(x,y,z)=esz(x,y,z)-dts*(DX(e1y)-DY(e1x))
#endif
  enddo; enddo; enddo
!     invert Helmholtz equation to obtain B2
#ifdef ELM
  call helmholtz(b2x,esx,real(eps,drk),fixed_cycles) 
  call helmholtz(b2y,esy,real(eps,drk),fixed_cycles)
  call helmholtz(b2z,esz,real(eps,drk),fixed_cycles)
#else
  b2x=esx; b2y=esy; b2z=esz
#endif

    do z=1,nz; do y=1,ny; do x=1,nx
      jt1x(x,y,z) = -DELTA(b1x)
      jt1y(x,y,z) = -DELTA(b1y)
      jt1z(x,y,z) = -DELTA(b1z)
    enddo; enddo; enddo
    call bound_j(jt1x,jt1y,jt1z)
#ifdef EXPANSION
    do z=1,nz; do y=1,ny; do x=1,nx
      b2x(x,y,z)=b2x(x,y,z)+dts*(nu*DELTA(jt1x)/n_0-eta*jt1x(x,y,z)-     2*daoa*b1x(x,y,z))
      b2y(x,y,z)=b2y(x,y,z)+dts*(nu*DELTA(jt1y)/n_0-eta*jt1y(x,y,z)-(daoa+dwow)*b1y(x,y,z))
      b2z(x,y,z)=b2z(x,y,z)+dts*(nu*DELTA(jt1z)/n_0-eta*jt1z(x,y,z)-(daoa+dwow)*b1z(x,y,z))
    enddo; enddo; enddo
#else
    do z=1,nz; do y=1,ny; do x=1,nx
      b2x(x,y,z)=b2x(x,y,z)+dts*(nu*DELTA(jt1x)/n_0-eta*jt1x(x,y,z))
      b2y(x,y,z)=b2y(x,y,z)+dts*(nu*DELTA(jt1y)/n_0-eta*jt1y(x,y,z))
      b2z(x,y,z)=b2z(x,y,z)+dts*(nu*DELTA(jt1z)/n_0-eta*jt1z(x,y,z))
    enddo; enddo; enddo
#endif

#ifdef NO_ELEC_PHYS
#else
#ifndef isothermal
#ifndef barotropic
! Calculate electron velocity
  call bound_j(jx,jy,jz)
  call bound_j(jtx,jty,jtz)
  call bound_rho(rho)
! corrected the bug that dts below had a positive sign
! when it should have been negative. Thanks to kittipat
! 		Tulasi: 10/22/2007
  do z=1,nz; do y=1,ny; do x=1,nx
    vex(x,y,z)=(jx(x,y,z)-jtx(x,y,z))/rho(x,y,z) 
    vey(x,y,z)=(jy(x,y,z)-jty(x,y,z))/rho(x,y,z) 
    vez(x,y,z)=(jz(x,y,z)-jtz(x,y,z))/rho(x,y,z)
  enddo; enddo; enddo
  call bound_j(vex,vey,vez)
    te=pe2/rho
!     advance pe2
  do z=1,nz; do y=1,ny; do x=1,nx
    jt1x(x,y,z)=vex(x,y,z)*pe1(x,y,z)*rho(x,y,z)**(1-gam)
    jt1y(x,y,z)=vey(x,y,z)*pe1(x,y,z)*rho(x,y,z)**(1-gam)
    jt1z(x,y,z)=vez(x,y,z)*pe1(x,y,z)*rho(x,y,z)**(1-gam)
    pe2(x,y,z)=pe2(x,y,z)*rho(x,y,z)**(1-gam)
  enddo; enddo; enddo
  do z=1,nz; do y=1,ny; do x=1,nx
    jt2x(x,y,z)=DELTA(b1x)
    jt2y(x,y,z)=DELTA(b1y)
    jt2z(x,y,z)=DELTA(b1z)
  enddo; enddo; enddo
  call boundary(0,0,jt1x)
  call boundary(0,0,jt1y)
  call boundary(0,0,jt1z)
! added the $$-\nabla\dot q + \alpha \nu |\nabla\times J|^2$$ to the RHS
#ifdef EXPANSION
  do z=1,nz; do y=1,ny; do x=1,nx
    pe2(x,y,z)=pe2(x,y,z)-dts*(DX(jt1x)+DY(jt1y)+DZ(jt1z) - &
      (2*daoa+dwow)*pe2(x,y,z) -&
      (gam-1)*(chi*DELTA(te)- &
      nu*((jt2x(x,y,z))**2 + (jt2y(x,y,z))**2 +(jt2z(x,y,z))**2)/n_0&
    )/rho(x,y,z)**(gam-1))
  enddo; enddo; enddo
#else
  do z=1,nz; do y=1,ny; do x=1,nx
    pe2(x,y,z)=pe2(x,y,z)-dts*(DX(jt1x)+DY(jt1y)+DZ(jt1z) - &
      (gam-1)*(chi*DELTA(te)- &
      nu*((jt2x(x,y,z))**2 + (jt2y(x,y,z))**2 +(jt2z(x,y,z))**2)/n_0&
    )/rho(x,y,z)**(gam-1))
  enddo; enddo; enddo
#endif

!get P back from the P/n^(g-1)
  do z=1,nz; do y=1,ny; do x=1,nx
    pe2(x,y,z)=pe2(x,y,z)*rho(x,y,z)**(gam-1)
  enddo; enddo; enddo
  call bound_pe(pe2)
#endif
#endif
#endif
  call bound_b(b2x,b2y,b2z)
!Calculate J_tot
  do z=1,nz; do y=1,ny; do x=1,nx
    jtx(x,y,z)=DY(b2z)-DZ(b2y) 
    jty(x,y,z)=DZ(b2x)-DX(b2z) 
    jtz(x,y,z)=DX(b2y)-DY(b2x)
  enddo; enddo; enddo

! second half step: b1,pe1 t --> t+dts
!     calculate E = ((curl b)x(1-d_e^2 delt^2)B - (ji x B) - delt pe)/n
  do z=1,nz; do y=1,ny; do x=1,nx
    esx(x,y,z)=HELMHOLTZ(b2x) 
    esy(x,y,z)=HELMHOLTZ(b2y) 
    esz(x,y,z)=HELMHOLTZ(b2z)
    e1x(x,y,z)=(VECX(jt,es)-VECX(j,b2))/rho(x,y,z) 
    e1y(x,y,z)=(VECY(jt,es)-VECY(j,b2))/rho(x,y,z) 
    e1z(x,y,z)=(VECZ(jt,es)-VECZ(j,b2))/rho(x,y,z) 
!   e1x(x,y,z)=(V_ECX(jt,es)-V_ECX(j,b2)-aleph*D_X(pe2))/rho(x,y,z) 
!   e1y(x,y,z)=(V_ECY(jt,es)-V_ECY(j,b2)-aleph*D_Y(pe2))/rho(x,y,z) 
!   e1z(x,y,z)=(V_ECZ(jt,es)-V_ECZ(j,b2)-aleph*D_Z(pe2))/rho(x,y,z) 
  enddo; enddo; enddo
!     calculate B1' from B1, advance B1'
  call bound_e(e1x,e1y,e1z)
  do z=1,nz; do y=1,ny; do x=1,nx
#ifdef DRIVE_B
    esx(x,y,z)=HELMHOLTZ(b1x); esx(x,y,z)=esx(x,y,z)-dts*(DY(e1z)-DZ(e1y)-cos(FOUROMEGA*t)*dr_bx(x,y,z))
    esy(x,y,z)=HELMHOLTZ(b1y); esy(x,y,z)=esy(x,y,z)-dts*(DZ(e1x)-DX(e1z)-cos(FOUROMEGA*t)*dr_by(x,y,z))
    esz(x,y,z)=HELMHOLTZ(b1z); esz(x,y,z)=esz(x,y,z)-dts*(DX(e1y)-DY(e1x)-cos(FOUROMEGA*t)*dr_bz(x,y,z))

#else
    esx(x,y,z)=HELMHOLTZ(b1x); esx(x,y,z)=esx(x,y,z)-dts*(DY(e1z)-DZ(e1y))
    esy(x,y,z)=HELMHOLTZ(b1y); esy(x,y,z)=esy(x,y,z)-dts*(DZ(e1x)-DX(e1z))
    esz(x,y,z)=HELMHOLTZ(b1z); esz(x,y,z)=esz(x,y,z)-dts*(DX(e1y)-DY(e1x))
#endif
  enddo; enddo; enddo
! Calculate \nabla^2(b2x) before overwriting B2. Save it to jt1 
    do z=1,nz; do y=1,ny; do x=1,nx
      jt1x(x,y,z) = -DELTA(b2x)
      jt1y(x,y,z) = -DELTA(b2y)
      jt1z(x,y,z) = -DELTA(b2z)
    enddo; enddo; enddo
!     save B1 to connect leapfrog grids
  do z=1,nz; do y=1,ny; do x=1,nx
    b2x(x,y,z)=b1x(x,y,z); b2y(x,y,z)=b1y(x,y,z); b2z(x,y,z)=b1z(x,y,z)
  enddo; enddo; enddo
!     invert Helmholtz equation to obtain B1
#ifdef ELM
  call helmholtz(b1x,esx,real(eps,drk),fixed_cycles) 
  call helmholtz(b1y,esy,real(eps,drk),fixed_cycles)
  call helmholtz(b1z,esz,real(eps,drk),fixed_cycles)
#else
  b1x=esx; b1y=esy; b1z=esz
#endif

    call bound_j(jt1x,jt1y,jt1z)
#ifdef EXPANSION
    do z=1,nz; do y=1,ny; do x=1,nx
      b1x(x,y,z)=b1x(x,y,z)+dts*(nu*DELTA(jt1x)/n_0-eta*jt1x(x,y,z)-     2*daoa*b2x(x,y,z))
      b1y(x,y,z)=b1y(x,y,z)+dts*(nu*DELTA(jt1y)/n_0-eta*jt1y(x,y,z)-(daoa+dwow)*b2y(x,y,z))
      b1z(x,y,z)=b1z(x,y,z)+dts*(nu*DELTA(jt1z)/n_0-eta*jt1z(x,y,z)-(daoa+dwow)*b2z(x,y,z))
    enddo; enddo; enddo
#else
    do z=1,nz; do y=1,ny; do x=1,nx
      b1x(x,y,z)=b1x(x,y,z)+dts*(nu*DELTA(jt1x)/n_0-eta*jt1x(x,y,z))
      b1y(x,y,z)=b1y(x,y,z)+dts*(nu*DELTA(jt1y)/n_0-eta*jt1y(x,y,z))
      b1z(x,y,z)=b1z(x,y,z)+dts*(nu*DELTA(jt1z)/n_0-eta*jt1z(x,y,z))
    enddo; enddo; enddo
#endif


#ifdef NO_ELEC_PHYS
#else
#ifndef isothermal
#ifndef barotropic
! Calculate electron velocity
  call bound_j(jx,jy,jz)
  call bound_j(jtx,jty,jtz)
  call bound_rho(rho)
! corrected the bug that dts below had a positive sign
! when it should have been negative. Thanks to kittipat
! 		Tulasi: 10/22/2007
  do z=1,nz; do y=1,ny; do x=1,nx
    vex(x,y,z)=(jx(x,y,z)-jtx(x,y,z))/rho(x,y,z) 
    vey(x,y,z)=(jy(x,y,z)-jty(x,y,z))/rho(x,y,z) 
    vez(x,y,z)=(jz(x,y,z)-jtz(x,y,z))/rho(x,y,z)
  enddo; enddo; enddo
  call bound_j(vex,vey,vez)
    te=pe1/rho
!     save pe1 to connect grids, advance pe1 
  do z=1,nz; do y=1,ny; do x=1,nx
    esx(x,y,z) = pe1(x,y,z)
    jt1x(x,y,z)=vex(x,y,z)*pe2(x,y,z)*rho(x,y,z)**(1-gam)
    jt1y(x,y,z)=vey(x,y,z)*pe2(x,y,z)*rho(x,y,z)**(1-gam)
    jt1z(x,y,z)=vez(x,y,z)*pe2(x,y,z)*rho(x,y,z)**(1-gam)
    pe1(x,y,z)=pe1(x,y,z)*rho(x,y,z)**(1-gam)
  enddo; enddo; enddo
  do z=1,nz; do y=1,ny; do x=1,nx
    jt2x(x,y,z)=DELTA(b2x)
    jt2y(x,y,z)=DELTA(b2y)
    jt2z(x,y,z)=DELTA(b2z)
  enddo; enddo; enddo
  call boundary(0,0,jt1x)
  call boundary(0,0,jt1y)
  call boundary(0,0,jt1z)
#ifdef EXPANSION
  do z=1,nz; do y=1,ny; do x=1,nx
    pe1(x,y,z)=pe1(x,y,z)-dts*(DX(jt1x)+DY(jt1y)+DZ(jt1z) - &
      (2*daoa+dwow)*pe1(x,y,z) - &
      (gam-1)*(chi*DELTA(te) - &
      nu*((jt2x(x,y,z))**2 + (jt2y(x,y,z))**2 +(jt2z(x,y,z))**2)/n_0&
    )/rho(x,y,z)**(gam-1))
  enddo; enddo; enddo
#else
  do z=1,nz; do y=1,ny; do x=1,nx
    pe1(x,y,z)=pe1(x,y,z)-dts*(DX(jt1x)+DY(jt1y)+DZ(jt1z) - &
      (gam-1)*(chi*DELTA(te) - &
      nu*((jt2x(x,y,z))**2 + (jt2y(x,y,z))**2 +(jt2z(x,y,z))**2)/n_0&
    )/rho(x,y,z)**(gam-1))
  enddo; enddo; enddo
#endif

! get P back from the P/n^(g-1)
  do z=1,nz; do y=1,ny; do x=1,nx
    pe1(x,y,z)=pe1(x,y,z)*rho(x,y,z)**(gam-1)
  enddo; enddo; enddo
#endif
#endif
#endif
! connect leapfrog grids
  do z=1,nz; do y=1,ny; do x=1,nx
    b2x(x,y,z) = (b2x(x,y,z)+b1x(x,y,z))/2
    b2y(x,y,z) = (b2y(x,y,z)+b1y(x,y,z))/2
    b2z(x,y,z) = (b2z(x,y,z)+b1z(x,y,z))/2
#ifdef NO_ELEC_PHYS
#else
#ifndef isothermal
#ifndef barotropic
    pe2(x,y,z) = (esx(x,y,z)+pe1(x,y,z))/2
#endif
#endif
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
#ifdef EXPANSION
  Use expansion_params
#endif
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
#ifdef EXPANSION
  Use expansion_params
#endif
  Implicit None

  Integer x,y,z
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), & 
                               dz=lz/(nz*pez), dts=dt/(2*substeps)
! Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: tpe, tro

  call start_timer()
  call bound_b(b1x,b1y,b1z)
  call bound_pe(pe1)
! call smooth(pe1,tpe)
! call smooth(rho,tro)
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
