!******************************************************************************
!          3D particle code: fluid equations for two fluid version
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
#ifndef D_pi
#define D_pi 0.
#endif
#ifndef D_j
#define D_j 0.
#endif
#ifndef D_n
#define D_n 0.
#endif
#ifndef n_base
#define n_base 0.
#endif

#define X x,y,z
#define DOALL do z=1,nz; do y=1,ny; do x=1,nx
#define DONE enddo; enddo; enddo

define(_,`
  ')
! XYZ: generate argument three times with $1 substituted by x,y,z:
!   e.g.: XYZ(jt$1(X)=CURL$1(b1)) -> jtx(X)=CURLx(b1); jty(X)=CURLy(b1); ...
! SUM: like XYZ but add the three results instead of creating commands
!   e.g.: SUM(j1$1(X)*D$1(tpy)) -> j1x(X)*Dx(tpy)+j1y(X)*Dy(tpy)+j1z(X)*Dz(tpy)
define(XYZ,`define(AA,`$1') AA(x) _ AA(y) _ AA(z) undefine(`AA')')
define(SUM,`define(AA,`$1') & _ (AA(x)+& _ AA(y)+& _ AA(z)) undefine(`AA')')
define(DELTA,`(&
  ($1(x+1,y,z)-2*$1(x,y,z)+$1(x-1,y,z))/dx**2+&
  ($1(x,y+1,z)-2*$1(x,y,z)+$1(x,y-1,z))/dy**2+&
  ($1(x,y,z+1)-2*$1(x,y,z)+$1(x,y,z-1))/dz**2& _ )')
define(HELMHOLTZ,`(& _ $1(x,y,z)-d_e2*DELTA($1))')
define(Dx,`(& _ ($1(x+1,y,z)-$1(x-1,y,z))/(2*dx)& _ )')
define(Dy,`(& _ ($1(x,y+1,z)-$1(x,y-1,z))/(2*dy)& _ )')
define(Dz,`(& _ ($1(x,y,z+1)-$1(x,y,z-1))/(2*dz)& _ )')
define(VECx,`(& _ $1y(x,y,z)*$2z(x,y,z)-$1z(x,y,z)*$2y(x,y,z)& _ )')
define(VECy,`(& _ $1z(x,y,z)*$2x(x,y,z)-$1x(x,y,z)*$2z(x,y,z)& _ )')
define(VECz,`(& _ $1x(x,y,z)*$2y(x,y,z)-$1y(x,y,z)*$2x(x,y,z)& _ )')
define(DIV,`(Dx($1x)+Dy($1y)+Dz($1z))')
define(CURLx,`(Dy($1z)-Dz($1y))')
define(CURLy,`(Dz($1x)-Dx($1z))')
define(CURLz,`(Dx($1y)-Dy($1x))')

!------------------------------------------------------------------------------
!                      step all fluid quantities
!------------------------------------------------------------------------------

! advance variables according to two-fluid equations (leapfrog algorithm)
! b1, j1, pi1, pe1, n1: values at time t (on entry); at time t+dt (on exit)
! b2, j2, pi2, pe2, n2: values at time t-dt/2 (on entry)
!                 average of values at time t and time t+dt/2 (on exit)

Subroutine steptwofluid(fixed_cycles)
  Use pe_env
  Use partfield
  Use timer
  Implicit None

  Integer x,y,z,i,dampcycles
  Logical :: fixed_cycles
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez)
!     jt: temporary array to store the total current j = curl B
!     bp: temporary array to store B-d_e^2 delt^2 B
!     et: temporary array to store the electric field
!     tp: general purpose temporary array
!     bs,js,pes,pis,ns: save arrays to connect leapfrog-grids
  Real(kind=drk), Save, Dimension(0:nx+1,0:ny+1,0:nz+1) :: &
        bpx,bpy,bpz,jtx,jty,jtz,etx,ety,etz,tpx,tpy,tpz, &
        bsx,bsy,bsz,jsx,jsy,jsz,pes,pis,ns

  call start_timer()
  DOALL; n1(X)=n1(X)+n_base; n2(X)=n2(X)+n_base; DONE
! first half step: b2,j2,pi2,pe2,n2 t-dt/2 --> t+dt/2
!     calculate E = ((curl b)x(1-d_e^2 delt^2)B - (ji x B) - delt pe)/n
  call bound_b(b1x,b1y,b1z); call bound_pe(pe1); call bound_pi(pi1)
  DOALL
    XYZ(bp$1(X)=HELMHOLTZ(b1$1)) ; XYZ(jt$1(X)=CURL$1(b1))
    XYZ(et$1(X)=(VEC$1(jt,bp)-VEC$1(j1,b1)-D$1(pe1))/n1(X))
  DONE
!     calculate B2' from B2, advance B2', invert Helmholtz equ. to obtain B2
  call bound_e(etx,ety,etz); call bound_b(b2x,b2y,b2z)
  DOALL; XYZ(bp$1(X)=HELMHOLTZ(b2$1)-dt*CURL$1(et)); DONE
  XYZ(call helmholtz(b2$1,bp$1,real(eps,drk),fixed_cycles))
!     advance ion current j2x,j2y,j2z
  call bound_j(j1x,j1y,j1z); call bound_n(n1)
  XYZ(tpx=j1x/n1*j1$1; tpy=j1y/n1*j1$1; tpz=j1z/n1*j1$1
      DOALL; j2$1(X)=j2$1(X)-dt*(DIV(tp)+D$1(pe1)+D$1(pi1)-VEC$1(jt,b1)); DONE)
!     advance n2, pe2, pi2 in conservative form, second order viscosities
#ifndef isothermal
  tpx=pe1*n1**(-gam_e); tpy=pi1*n1**(-gam_i)
  DOALL
    pe2(X)=pe2(X)*n2(X)**(-gam_e)-dt*1/n1(X)*SUM((j1$1(X)-jt$1(X))*D$1(tpx))
    pi2(X)=pi2(X)*n2(X)**(-gam_i)-dt*1/n1(X)*SUM(j1$1(X)*D$1(tpy))
    n2(X)=abs(n2(X)-dt*DIV(j1)-n_base)+n_base
    pe2(X)=pe2(X)*n2(X)**(gam_e)
    pi2(X)=pi2(X)*n2(X)**(gam_i)
  DONE
#else
  DOALL
    n2(X)=abs(n2(X)-dt*DIV(j1)-n_base)+n_base
    pe2(X)=T_e*n2(X)
    pi2(X)=T_i*n2(X)
  DONE
#endif
! fourth order hyperviscosity for B, j, pe, pi, n
  dampcycles=int(max(D_b,D_j,D_n,D_pe,D_pi)*dt*(3.142/min(dx,dy,dz))**4*20.+1.)
  do i=1,dampcycles
    call bound_b(b2x,b2y,b2z); call bound_j(j2x,j2y,j2z)
    DOALL; XYZ(tp$1(X)=DELTA(b2$1)); XYZ(jt$1(X)=DELTA(j2$1)); DONE
    call bound_b(tpx,tpy,tpz); call bound_j(jtx,jty,jtz)
    DOALL 
      XYZ(b2$1(X)=b2$1(X)-D_b*dt/dampcycles*DELTA(tp$1))
      XYZ(j2$1(X)=j2$1(X)-D_j*dt/dampcycles*DELTA(jt$1))
    DONE
    call bound_n(n2); call bound_pe(pe2); call bound_pi(pi2)
    DOALL; tpx(X)=DELTA(n2); tpy(X)=DELTA(pe2); tpz(X)=DELTA(pi2); DONE
    call bound_n(tpx); call bound_pe(tpy); call bound_pi(tpz)
    DOALL
      n2(X)=abs(n2(X)-D_n*dt/dampcycles*DELTA(tpx)-n_base)+n_base
#ifndef isothermal
      pe2(X)=pe2(X)-D_pe*dt/dampcycles*DELTA(tpy)
      pi2(X)=pi2(X)-D_pi*dt/dampcycles*DELTA(tpz)
#endif
    DONE
  enddo
! second half step: b1,j1,pi1,pe1,n1 t --> t+dt
!     save old values to connect leapfrog-grids
  DOALL
    XYZ(bs$1(X)=b1$1(X)); XYZ(js$1(X)=j1$1(X))
    ns(X) = n1(X); pes(X)=pe1(X); pis(X)=pi1(X) 
  DONE
!     calculate E = ((curl b)x(1-d_e^2 delt^2)B - (ji x B) - delt pe)/n
  call bound_b(b2x,b2y,b2z); call bound_pe(pe2); call bound_pi(pi2)
  DOALL
    XYZ(bp$1(X)=HELMHOLTZ(b2$1)) ; XYZ(jt$1(X)=CURL$1(b2))
    XYZ(et$1(X)=(VEC$1(jt,bp)-VEC$1(j2,b2)-D$1(pe2))/n2(X))
  DONE
!     calculate B1' from B1, advance B1', invert Helmholtz equ. to obtain B1
  call bound_e(etx,ety,etz); call bound_b(b1x,b1y,b1z)
  DOALL; XYZ(bp$1(X)=HELMHOLTZ(b1$1)-dt*CURL$1(et)); DONE
  XYZ(call helmholtz(b1$1,bp$1,real(eps,drk),fixed_cycles))
!     advance ion current j1x,j1y,j1z
  call bound_j(j2x,j2y,j2z); call bound_n(n2)
  XYZ(tpx=j2x/n2*j2$1; tpy=j2y/n2*j2$1; tpz=j2z/n2*j2$1
      DOALL; j1$1(X)=j1$1(X)-dt*(DIV(tp)+D$1(pe2)+D$1(pi2)-VEC$1(jt,b2)); DONE)
!     advance n1, pe1, pi1 in conservative form, second order viscosities
#ifndef isothermal
  tpx=pe2*n2**(-gam_e); tpy=pi2*n2**(-gam_i)
  DOALL
    pe1(X)=pe1(X)*n1(X)**(-gam_e)-dt*1/n2(X)*SUM((j2$1(X)-jt$1(X))*D$1(tpx))
    pi1(X)=pi1(X)*n1(X)**(-gam_i)-dt*1/n2(X)*SUM(j2$1(X)*D$1(tpy))
    n1(X)=abs(n1(X)-dt*DIV(j2)-n_base)+n_base
    pe1(X)=pe1(X)*n1(X)**(gam_e)
    pi1(X)=pi1(X)*n1(X)**(gam_i)
  DONE
#else
  DOALL
    n1(X)=abs(n1(X)-dt*DIV(j2)-n_base)+n_base
    pe1(X)=T_e*n1(X)
    pi1(X)=T_i*n1(X)
  DONE
#endif
! fourth order hyperviscosity for B, j, pe, pi, n
  dampcycles=int(max(D_b,D_j,D_n,D_pe,D_pi)*dt*(3.142/min(dx,dy,dz))**4*20.+1.)
  do i=1,dampcycles
    call bound_b(b1x,b1y,b1z); call bound_j(j1x,j1y,j1z)
    DOALL; XYZ(tp$1(X)=DELTA(b1$1)); XYZ(jt$1(X)=DELTA(j1$1)); DONE
    call bound_b(tpx,tpy,tpz); call bound_j(jtx,jty,jtz)
    DOALL 
      XYZ(b1$1(X)=b1$1(X)-D_b*dt/dampcycles*DELTA(tp$1))
      XYZ(j1$1(X)=j1$1(X)-D_j*dt/dampcycles*DELTA(jt$1))
    DONE
    call bound_n(n1); call bound_pe(pe1); call bound_pi(pi1)
    DOALL; tpx(X)=DELTA(n1); tpy(X)=DELTA(pe1); tpz(X)=DELTA(pi1); DONE
    call bound_n(tpx); call bound_pe(tpy); call bound_pi(tpz)
    DOALL
      n1(X)=abs(n1(X)-D_n*dt/dampcycles*DELTA(tpx)-n_base)+n_base
#ifndef isothermal
      pe1(X)=pe1(X)-D_pe*dt/dampcycles*DELTA(tpy)
      pi1(X)=pi1(X)-D_pi*dt/dampcycles*DELTA(tpz)
#endif
    DONE
  enddo
!     connect leapfrog-grids
  DOALL
    XYZ(b2$1(X)=(bs$1(X)+b1$1(X))/2) _ XYZ(j2$1(X)=(js$1(X)+j1$1(X))/2)
    n2(X)=(ns(X)+n1(X))/2 _ pe2(X)=(pes(X)+pe1(X))/2 
    pi2(X)=(pis(X)+pi1(X))/2
  DONE
  DOALL; n1(X)=n1(X)-n_base; n2(X)=n2(X)-n_base; DONE
  call stop_timer(t_step2fluid,n_step2fluid)

End Subroutine steptwofluid

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
  Real(kind=drk), Save, Dimension(0:nx+1,0:ny+1,0:nz+1) :: tmp1, tmp2

  tmp1(:,:,:) = 0
  call bound_b(bx,by,bz)
  DOALL; tmp2(X)=-DIV(b); DONE
  call poisson(tmp1,tmp2,real(eps,drk))
  call bound_phi(0,0,tmp1)
  DOALL; XYZ(b$1(X)=b$1(X)+D$1(tmp1)); DONE
End Subroutine divergence
