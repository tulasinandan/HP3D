!******************************************************************************
!          3D particle code: Miscelleneous Subroutines written by Tulasi 
!                          Tulasi Nandan Parashar, 2006 
!
!                       LATEST CHANGE: November 13, 2006
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


!-----------------------------------------------------------------
! Subroutine to calculate the fluid elctron current from the ion
! current and the magnetic field.
!					Tulasi Nandan Parashar
!					11/13/2006
!-----------------------------------------------------------------
Subroutine calc_je(bxx,byy,bzz,jix,jiy,jiz,jelecx,jelecy,jelecz,jtox,jtoy,jtoz)
    use pe_env
    use partfield
    use timer
#ifdef EXPANSION
    Use expansion_params
#endif

    Implicit None
    Integer :: x,y,z
    Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez)
    Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: bxx,byy,bzz,jix,jiy,jiz,jelecx,jelecy,jelecz
    Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: jtox,jtoy,jtoz
    Logical :: reset

    call start_timer()
    
  call bound_b(bxx,byy,bzz)
    do z=1,nz; do y=1,ny; do x=1,nx 
      jtox(x,y,z)=DY(bzz)-DZ(byy)
      jtoy(x,y,z)=DZ(bxx)-DX(bzz)
      jtoz(x,y,z)=DX(byy)-DY(bxx)

      jelecx(x,y,z)=jix(x,y,z)-jtox(x,y,z) 
      jelecy(x,y,z)=jiy(x,y,z)-jtoy(x,y,z)
      jelecz(x,y,z)=jiz(x,y,z)-jtoz(x,y,z)
    enddo; enddo; enddo

    call bound_rho_j(jelecx); call bound_rho_j(jelecy); call bound_rho_j(jelecz) 
    call bound_rho_j(jtox); call bound_rho_j(jtoy); call bound_rho_j(jtoz) 
    call stop_timer(t_calc_je, n_calc_je)
End Subroutine calc_je



!-----------------------------------------------------------------
! Subroutine to calculate the vorticity of a field.
!					Tulasi Nandan Parashar
!					06/29/2007
!-----------------------------------------------------------------
!Subroutine vorticity(fieldx,fieldy,fieldz,vortix,vortiy,vortiz)
!    use pe_env
!    use partfield
!    use timer
!
!    Implicit None
!    Integer :: x,y,z
!    Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez)
!    Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: fieldx,fieldy,fieldz,vortix,vortiy,vortiz
!    Logical :: reset
!
!    call start_timer()
!    
!    do z=1,nz; do y=1,ny; do x=1,nx 
!      vortix(x,y,z)=DDDY(fieldz)-DDDZ(fieldy)
!      vortiy(x,y,z)=DDDZ(fieldx)-DDDX(fieldz)
!      vortiz(x,y,z)=DDDX(fieldy)-DDDY(fieldx)
!    enddo; enddo; enddo
!
!    call bound_rho_j(vortix); call bound_rho_j(vortiy); call bound_rho_j(vortiz) 
!    call stop_timer(t_vorticity, n_vorticity)
!End Subroutine vorticity


 Subroutine set_drive_fields()
   use pe_env
   use partfield
   use timer
   Implicit None

   Real(kind=8), Dimension(1:nx*pex,1:ny*pey,1:nz*pez) :: tmp
   integer :: nnnn
   Integer :: i,j,k
   inquire(iolength=nnnn) tmp
#ifdef DRIVE_B
    if (myproc==0) print *, 'Reached inside set drive fields'
    i=my_pex*nx; j=my_pey*ny; k=my_pez*nz
    print *, myproc, i, j, k
    open(unit=8,status='old',file='Fbx.dat',form='unformatted',access='direct',recl=nnnn)
    read(8,rec=1) tmp
    close(8)
    dr_bx(:,:,:) = alpha_drb*tmp(i+1:i+nx,j+1:j+ny,k+1:k+nz)
    open(unit=8,status='old',file='Fby.dat',form='unformatted',access='direct',recl=nnnn)    
    read(8,rec=1) tmp
    close(8)
    dr_by(:,:,:) = alpha_drb*tmp(i+1:i+nx,j+1:j+ny,k+1:k+nz)
    open(unit=8,status='old',file='Fbz.dat',form='unformatted',access='direct',recl=nnnn)    
    read(8,rec=1) tmp
    close(8)
    dr_bz(:,:,:) = alpha_drb*tmp(i+1:i+nx,j+1:j+ny,k+1:k+nz)
    if (myproc==0) print *, 'read drive fields'
#else
#endif
 end Subroutine set_drive_fields
!-----------------------------------------------------------------------------------------------
! Calculate different quantities for analysis
!                    T. N. Parashar - 03/11/2009
!-----------------------------------------------------------------------------------------------

!Subroutine analyze(nstp1,nstp2)

!   use pe_env
!   use partfield
!   use timer

!     Implicit None
!     Real(kind=8)                    :: nstp1, nstp2
!     Real(kind=8), Dimension(0:nx+1,0:ny+1,0:nz+1)  :: vix, viy, viz, energyj, vdjcb
!     Real(kind=8), Dimension(0:nx+1,0:ny+1,0:nz+1)  :: vdgp,venstrz, energyi, energyb
!     Real(kind=8), Dimension(nstp2)                 :: meantx, meanty, meantz, gpterm, dekdt, dbsqdt
!     Real(kind=8), Dimension(nstp2)                 :: djsqdt, vcjterm, meanoz, meanjzsq, avgenergybulk
!     Real(kind=8), Dimension(nstp2)                 :: avgenergyb, avgenergyj

!   call start_timer()

!        vix=jx/rho; viy=jy/rho; viz=jz/rho 

!        pxx=smooth_t(pxx,5); pyy=smooth_t(pyy,5); pzz=smooth_t(pzz,5)
!        jtotx=smooth_t(jtotx,40); jtoty=smooth_t(jtoty,40); jtotz=smooth_t(jtotz,40)
!        vix=smooth_t(vix,40); viy=smooth_t(viy,40); viz=smooth_t(viz,40);n=smooth_t(n,5)
!        
!        vdgp=vix*deriv((pxx+pyy+pzz)/3.,xx,1)+viy*deriv((pxx+pyy+pzz)/3.,yy,2)
!        energyi=(vix**2+viy**2+viz**2)*n/2; energyb=(bx1*bx1+by*by)/2.
!        energyj=0.5*d_e2*(jtotx**2+jtoty**2+jtotz**2)

!        vdjcb=vix*(jtoty*bz1-jtotz*by1)+viy*(jtotz*bx1-jtotx*bz1)+viz*(jtotx*by1-jtoty*bx1)
!        venstrz=(deriv(viy,xx,1)-deriv(vix,yy,2))**2

!        vcjterm(t)=mean(vdjcb); gpterm(t)=mean(vdgp); meanjzsq(t)=mean(jtotz**2)
!        avgenergybulk(t)=mean(energyi); avgenergyb(t)=mean(energyb); avgenergyj(t)=mean(energyj)
!        meantx(t)=mean(pxx/n); meanty(t)=mean(pyy/n); meantz(t)=mean(pzz/n); meanoz(t)=mean(venstrz)
!     if(nstp1 == nstp2) then
!     do t=2,Nsteps-1
!      dekdt(t)=(avgenergybulk(t+1)-avgenergybulk(t-1))/(2.*movieout)
!      dbsqdt(t)=(avgenergyb(t+1)-avgenergyb(t-1))/(2.*movieout)
!      djsqdt(t)=(avgenergyj(t+1)-avgenergyj(t-1))/(2.*movieout)
!     enddo
!     endif
!  CONTAINS
!     Function smooth_t(a,n)
!        Real(kind=8), Dimension(:,:,:) :: a
!        Integer            :: n
!        Integer                        :: i,j,k,ax,ay
!        Real(kind=8), Dimension(size(a,1),size(a,2),size(a,3)) :: smooth_t
!        Real(kind=8), Dimension(0:size(a,1)+1,0:size(a,2)+1,0:size(a,3)+1) :: tmp
!          ax=size(a,1); ay=size(a,2)
!          tmp=0.; tmp(1:ax,1:ay,:) = a(:,:,:)
!           tmp(0,1:ay,:) = a(ax,:,:); tmp(ax+1,1:ay,:) = a(1,:,:) 
!           tmp(1:ax,0,:) = a(:,ay,:); tmp(1:ax,ay+1,:) = a(:,1,:)
!        do k=1,n
!           do j=1,ay; do i=1,ax
!           tmp(i,j,:)=0.25*tmp(i,j,:)+0.125*(tmp(i-1,j,:)+tmp(i,j-1,:)+tmp(i+1,j,:)+tmp(i,j+1,:))&
!                    +0.0625*(tmp(i-1,j-1,:)+tmp(i-1,j+1,:)+tmp(i+1,j-1,:)+tmp(i+1,j+1,:))
!           enddo; enddo
!           tmp(0,:,:) = tmp(ax,:,:); tmp(ax+1,:,:) = tmp(1,:,:) 
!           tmp(:,0,:) = tmp(:,ay,:); tmp(:,ay+1,:) = tmp(:,1,:)
!        enddo
!        smooth_t=tmp(0:ax,0:ay,:)
!     end function smooth_t

!     Function deriv(a,d,n)
!        Implicit None
!        Real(kind=8), Dimension(:,:,:)     :: a
!        Integer                              :: n
!        Real(kind=8), Dimension(size(a,n))   :: d
!        Real(kind=8), Dimension(0:size(a,1)+1,0:size(a,2)+1,0:size(a,3)+1) :: tmp
!        Integer                                          :: ax,ay
!        Real*8                                           :: delta
!        Real(kind=8), Dimension(size(a,1),size(a,2))     :: deriv
!          ax=size(a,1); ay=size(a,2)
!          tmp=0.; tmp(1:ax,1:ay,:) = a(:,:,:)
!          delta=d(2)-d(1)
!          if (n==1) then
!           tmp(0,1:ay,:) = a(ax,:,:); tmp(ax+1,1:ay,:) = a(1,:,:) 
!           tmp(1:ax,1:ay,:)=(tmp(2:ax+1,1:ay,:)-tmp(0:ax-1,1:ay,:))/(2*delta)
!          elseif (n==2) then
!           tmp(1:ax,0,:) = a(:,ay,:); tmp(1:ax,ay+1,:) = a(:,1,:)
!           tmp(1:ax,1:ay,:)=(tmp(1:ax,2:ay+1,:)-tmp(1:ax,0:ay-1,:))/(2*delta)
!          else
!           write(6,*) 'No direction specified for derivative. STOPPING'; STOP
!          endif
!          deriv = tmp(1:ax,1:ay,:)
!     end function deriv

!     Function mean(a)
!        Implicit None
!        Real(kind=8), Dimension(:,:,:) :: a
!        Real*8                       :: mean
!        mean=sum(a)/(size(a,1)*size(a,2)*size(a,3))
!     end function mean
!End Subroutine analyze
