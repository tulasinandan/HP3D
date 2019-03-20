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

#endif

!------------------------------------------------------------------------------'
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
  Real(kind=drk), Save, Dimension(0:nx+1,0:ny+1,0:nz+1) :: tnx,tny,tnz
#ifdef isothermal 
  Real(kind=drk), Save, Dimension(0:nx+1,0:ny+1,0:nz+1) :: kbx,kby,kbz
#elif barotropic
  Real(kind=drk), Save, Dimension(0:nx+1,0:ny+1,0:nz+1) :: kbx,kby,kbz
#else
  Real(kind=drk), Save, Dimension(0:nx+1,0:ny+1,0:nz+1) :: tnp,kbx,kby,kbz,kp,epe
#endif


  call start_timer()
#ifdef isothermal
    pe1=rho*T_e
#elif barotropic
    pe1=T_e*rho**gam
#else
#endif

#ifdef isothermal 
  esx=0.; esy=0.; esz=0.
#elif barotropic
  esx=0.; esy=0.; esz=0.
#else
  esx=0.; esy=0.; esz=0.; epe=0.
#endif

!!!##############################################
!!! CALCULATE K1=f(yn), ADD TO RHS
!!!##############################################
#ifdef isothermal 
  call calc_rk4_k(esx,esy,esz,pe1,kbx,kby,kbz)
#elif barotropic
  call calc_rk4_k(esx,esy,esz,pe1,kbx,kby,kbz)
#else
  call calc_rk4_k(esx,esy,esz,pe1,kbx,kby,kbz,kp)
#endif
  esx=esx+kbx/6.; tnx=b1x+dts*kbx/2.
  esy=esy+kby/6.; tny=b1y+dts*kby/2.
  esz=esz+kbz/6.; tnz=b1z+dts*kbz/2.
#ifdef isothermal 
#elif barotropic
#else
  epe=epe+kp/6. ; tnp=pe1+dts*kp/2.
#endif

!!!##############################################
!!! CALCULATE K2=f(yn+dts*K1/2.), ADD TO RHS
!!!##############################################
#ifdef isothermal 
  call calc_rk4_k(tnx,tny,tnz,pe1,kbx,kby,kbz)
#elif barotropic
  call calc_rk4_k(tnx,tny,tnz,pe1,kbx,kby,kbz)
#else
  call calc_rk4_k(tnx,tny,tnz,tnp,kbx,kby,kbz,kp)
#endif
  esx=esx+kbx/3.; tnx=b1x+dts*kbx/2.
  esy=esy+kby/3.; tny=b1y+dts*kby/2.
  esz=esz+kbz/3.; tnz=b1z+dts*kbz/2.
#ifdef isothermal 
#elif barotropic
#else
  epe=epe+kp/3. ; tnp=pe1+dts*kp/2.
#endif

!!!##############################################
!!! CALCULATE K3=f(yn+dts*K2/2.), ADD TO RHS
!!!##############################################
#ifdef isothermal 
  call calc_rk4_k(tnx,tny,tnz,pe1,kbx,kby,kbz)
#elif barotropic
  call calc_rk4_k(tnx,tny,tnz,pe1,kbx,kby,kbz)
#else
  call calc_rk4_k(tnx,tny,tnz,tnp,kbx,kby,kbz,kp)
#endif
  esx=esx+kbx/3.; tnx=b1x+dts*kbx 
  esy=esy+kby/3.; tny=b1y+dts*kby 
  esz=esz+kbz/3.; tnz=b1z+dts*kbz 
#ifdef isothermal 
#elif barotropic
#else
  epe=epe+kp/3. ; tnp=pe1+dts*kp 
#endif

!!!##############################################
!!! CALCULATE K4=f(yn+dts*K3), ADD TO RHS
!!!##############################################
#ifdef isothermal 
  call calc_rk4_k(tnx,tny,tnz,pe1,kbx,kby,kbz)
#elif barotropic
  call calc_rk4_k(tnx,tny,tnz,pe1,kbx,kby,kbz)
#else
  call calc_rk4_k(tnx,tny,tnz,tnp,kbx,kby,kbz,kp)
#endif
  esx=esx+kbx/6.
  esy=esy+kby/6.
  esz=esz+kbz/6.
#ifdef isothermal 
#elif barotropic
#else
  epe=epe+kp/6. 
#endif

!!!##############################################
!!! y(n+1) = y(n) + dts*(K1/6.+K2/3.+K3/3.+K4/6.) 
!!!##############################################
  b1x=b1x+dts*esx; b1y=b1y+dts*esy; b1z=b1z+dts*esz
#ifdef isothermal 
#elif barotropic
#else
  pe1=pe1+dts*epe
#endif

  call stop_timer(t_stephybrid,n_stephybrid)
End Subroutine stephybrid

!------------------------------------------------------------------------------
!                  CALCULATE THE K1, K2 etc. of RK4
!------------------------------------------------------------------------------
#ifdef isothermal 
Subroutine calc_rk4_k(tnx,tny,tnz,tnp,kbx,kby,kbz)
#elif barotropic
Subroutine calc_rk4_k(tnx,tny,tnz,tnp,kbx,kby,kbz)
#else
Subroutine calc_rk4_k(tnx,tny,tnz,tnp,kbx,kby,kbz,kp)
#endif
  Use pe_env
  Use partfield
#ifdef EXPANSION
  Use expansion_params
#endif
  Implicit None

  Integer :: x,y,z
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), & 
                               dz=lz/(nz*pez), dts=dt/(2*substeps)
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1), intent(in)  :: tnx, tny, tnz, tnp
#ifdef isothermal
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1), intent(out) :: kbx, kby, kbz
#elif barotropic
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1), intent(out) :: kbx, kby, kbz
#else
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1), intent(out) :: kbx, kby, kbz, kp
#endif

#ifdef isotropic
#elif barotropic
#else
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: te
#endif
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: jt1x,jt1y,jt1z

  call bound_pe(tnp)
  call bound_b(tnx,tny,tnz)
!Calculate J_tot
  do z=1,nz; do y=1,ny; do x=1,nx
    jt1x(x,y,z)=DY(tnz)-DZ(tny) 
    jt1y(x,y,z)=DZ(tnx)-DX(tnz) 
    jt1z(x,y,z)=DX(tny)-DY(tnx)
  enddo; enddo; enddo

!!!!! Calculate the electric field based on the following expression
!!!!!
!!!!! E' = J/n x es - Ji/n x b1 - grad(Pe)/n + eta* J + nu/n0 D_ELTA(J)
!!!!!

  do z=1,nz; do y=1,ny; do x=1,nx
    e1x(x,y,z)=(VECX(jt1,tn)-VECX(j,tn)-DX(tnp))/rho(x,y,z)+&
               eta*jt1x(x,y,z)+(nu/n_0)*DELTA(jt1x)
    e1y(x,y,z)=(VECY(jt1,tn)-VECY(j,tn)-DY(tnp))/rho(x,y,z)+&
               eta*jt1y(x,y,z)+(nu/n_0)*DELTA(jt1y)
    e1z(x,y,z)=(VECZ(jt1,tn)-VECZ(j,tn)-DZ(tnp))/rho(x,y,z)+&
               eta*jt1z(x,y,z)+(nu/n_0)*DELTA(jt1z)
  enddo; enddo; enddo

!!!!! Calculated kbx, kby, kbz etc. i.e. the RHS of dB/dt = f
! -NABLA X E'
  do z=1,nz; do y=1,ny; do x=1,nx
    kbx(x,y,z) = DZ(e1y)-DY(e1z)
    kby(x,y,z) = DX(e1z)-DZ(e1x)
    kbz(x,y,z) = DY(e1x)-DX(e1y)
  enddo; enddo; enddo
! ADD DRIVING TERM 
#ifdef DRIVE_B
    if (myproc == 0) write(6,*) 'cos(4wt)      ',cos(FOUROMEGA*t)
  do z=1,nz; do y=1,ny; do x=1,nx
    kbx(x,y,z) = kbx(x,y,z)-cos(FOUROMEGA*t)*dr_bx(x,y,z)
    kby(x,y,z) = kby(x,y,z)-cos(FOUROMEGA*t)*dr_by(x,y,z)
    kbz(x,y,z) = kbz(x,y,z)-cos(FOUROMEGA*t)*dr_bz(x,y,z)
  enddo; enddo; enddo
#endif
! ADD EXPANSION TERMS
#ifdef EXPANSION
  do z=1,nz; do y=1,ny; do x=1,nx
    kbx(x,y,z) = kbx(x,y,z)-     2*daoa*tnx(x,y,z)
    kby(x,y,z) = kby(x,y,z)-(daoa+dwow)*tny(x,y,z)
    kbz(x,y,z) = kbz(x,y,z)-(daoa+dwow)*tnz(x,y,z)
  enddo; enddo; enddo
#endif

#ifdef isothermal
#elif barotropic
#else
!!!############################################################################
!!! CALCULATE THE RHS OF dPe/dt = f(pe)
!!!
!!! where,
!!! f(pe)=gam*pe*grad.ve-ve.grad(pe)-n^gam*chi*grad^2(pe/n)-anisotropic pe+Q
!!! First implementation does not have anisotropic pe terms included.
!!!############################################################################
  te=tnp/rho
  do z=1,nz; do y=1,ny; do x=1,nx
    kp(x,y,z)=gam*tnp(x,y,z)*(DX(jt1x)-DX(jx)+DY(jt1y)-DY(jy)+DZ(jt1z)-DZ(jz))-&
    (jt1x(x,y,z)-jx(x,y,z))*DX(tnp)+(jt1y(x,y,z)-jy(x,y,z))*DY(tnp)+(jt1z(x,y,z)-&
    jz(x,y,z))*DZ(tnp)-chi*rho(x,y,z)**(gam)*DELTA(te)
  enddo; enddo; enddo
#endif

End Subroutine calc_rk4_k


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
