!******************************************************************************
! 2013/09/16 - Added the calculation of resistive energy loss to get a better
! handle on energy loss via resistivity. --> Tulasi 
!
!                   3D particle code: miscellaneous routines
!                         Andreas Zeiler, 1999-2001
!
!                        LATEST CHANGE: May 11, 2001
!
!******************************************************************************

#include "param"

#ifndef hybrid
#ifndef twofluid
#define fullparticle
#endif
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
!                                total energy
!------------------------------------------------------------------------------


Subroutine energy()
! calculate and print total energy denstity

  Use pe_env
  Use partfield
  Use timer
#ifdef EXPANSION
  Use expansion_params
#endif
  Implicit none

  Integer :: i,x,y,z
  Real(kind=drk) :: eem, ekin, eges, gamma, fact, ekinx, ekiny, ekinz
  Real(kind=drk) :: eelc, eprs, emag, eifl, eprsi, emagx, emagy, emagz 
  Real(kind=drk) :: eiflx, eifly, eiflz, eprsix, eprsiy, eprsiz, etaloss
  Real(kind=drk) :: elviscen
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez)
  Real(kind=drk) :: v_enstrophy_p, v_enstrophy_z, v_enstrophy_x, v_enstrophy_y
  Real(kind=drk) :: b_enstrophy_p, b_enstrophy_x, b_enstrophy_y, b_enstrophy_z
  Real(kind=drk), Save, Dimension(0:nx+1,0:ny+1,0:nz+1) :: vomegax,vomegay,vomegaz
  Real(kind=drk), Save, Dimension(0:nx+1,0:ny+1,0:nz+1) :: veloix,veloiy,veloiz
  Real(kind=drk), Save, Dimension(0:nx+1,0:ny+1,0:nz+1) :: bomegax,bomegay,bomegaz
  Real(kind=drk), Save, Dimension(0:nx+1,0:ny+1,0:nz+1) :: tmpx,tmpy,tmpz  
#ifdef twofluid
  Real(kind=drk), Save, Dimension(0:nx+1,0:ny+1,0:nz+1) :: jtx,jty,jtz
#endif

  call start_timer()

#ifdef fullparticle
! full particle version
!   electromagnetic field energy
  eem=0.
  eem = (sum(b1x(1:nx,1:ny,1:nz)**2)+sum(b1y(1:nx,1:ny,1:nz)**2)+&
        sum(b1z(1:nx,1:ny,1:nz)**2))/2+(sum(e1x(1:nx,1:ny,1:nz)**2)+&
        sum(e1y(1:nx,1:ny,1:nz)**2)+sum(e1z(1:nx,1:ny,1:nz)**2))/(2*c_2)
!   kinetic energy
  ekin=0.
#ifdef relativistic
  do i=1,np_i
    gamma=1/sqrt(1-(rvi(4,i)**2+rvi(5,i)**2+rvi(6,i)**2)/c_2)
    fact=gamma**2/(gamma+1)
    ekin=ekin+(rvi(4,i)**2+rvi(5,i)**2+rvi(6,i)**2)*fact/n_avg
  enddo
  do i=1,np_e
    gamma=1/sqrt(1-(rve(4,i)**2+rve(5,i)**2+rve(6,i)**2)/c_2)
    fact=gamma**2/(gamma+1)
    ekin=ekin+(rve(4,i)**2+rve(5,i)**2+rve(6,i)**2)*fact/n_avg*m_e
  enddo
#else
ekin = sum(rvi(4:6,1:np_i)**2)/(2*n_avg)
ekin = ekin + sum(rve(4:6,1:np_e)**2)/(2*n_avg)*m_e
#endif
!   total energy density
  eem=totalsum(eem)*(dx*dy*dz)/(lx*ly*lz) 
  ekin=totalsum(ekin)*(dx*dy*dz)/(lx*ly*lz) 
  eges=eem+ekin
#endif

#ifdef hybrid
! hybrid version
!   magnetic field energy, electron kinetic energy, electron thermal energy
!   1/2 B**2/2 + 1/2 d_e**2 j**2 + p_e/(gamma-1)
!   use esx, esy, esz as temporary arrays to store jx, jy, jz
  eem=0. ;  eelc=0.;  emag=0.;  eprs=0.; emagx=0.; emagy=0.; emagz=0.; elviscen=0.
    call bound_b(b1x,b1y,b1z)
    do z=1,nz; do y=1,ny; do x=1,nx
      tmpx(x,y,z) = DELTA(b1x)
      tmpy(x,y,z) = DELTA(b1y)
      tmpz(x,y,z) = DELTA(b1z)
    enddo; enddo; enddo
    do z=1,nz; do y=1,ny; do x=1,nx
      tmpx(x,y,z) = nu*tmpx(x,y,z)**2/n_0
      tmpy(x,y,z) = nu*tmpy(x,y,z)**2/n_0
      tmpz(x,y,z) = nu*tmpz(x,y,z)**2/n_0
    elviscen=tmpx(nx,ny,nz)+tmpy(nx,ny,nz)+tmpz(nx,ny,nz)
    enddo; enddo; enddo
!   elviscen=sum(tmpx(1:nx,1:ny,1:nz))+sum(tmpy(1:nx,1:ny,1:nz))+sum(tmpz(1:nx,1:ny,1:nz))
    elviscen=totalsum(elviscen)*(dx*dy*dz)/(lx*ly*lz)

  do z=1,nz; do y=1,ny; do x=1,nx
  if (rho(x,y,z) .lt. 1e-5) then
  write(6,*),'myproc',myproc,'xyz',x,y,z,rho(x,y,z)
  endif
  enddo; enddo; enddo

    call bound_rho(rho)
  do z=1,nz; do y=1,ny; do x=1,nx
    esx(x,y,z)=DY(b1z)-DZ(b1y)
    esy(x,y,z)=DZ(b1x)-DX(b1z)
    esz(x,y,z)=DX(b1y)-DY(b1x)
    etaloss = eta*(esx(x,y,z)**2+esy(x,y,z)**2+esz(x,y,z)**2)
    eem=eem+0.5*(b1x(x,y,z)**2+b1y(x,y,z)**2+b1z(x,y,z)**2)+ &
            0.5*d_e2*(esx(x,y,z)**2+esy(x,y,z)**2+esz(x,y,z)**2)+ &
            pe1(x,y,z)/(gam-1.)
    !
    eelc=eelc+0.5*d_e2*(esx(x,y,z)**2+esy(x,y,z)**2+esz(x,y,z)**2)
    eprs=eprs+pe1(x,y,z)/(gam-1.)
    emag=emag+0.5*(b1x(x,y,z)**2+b1y(x,y,z)**2+b1z(x,y,z)**2)
    emagx=emagx+0.5*b1x(x,y,z)**2
    emagy=emagy+0.5*b1y(x,y,z)**2
    emagz=emagz+0.5*b1z(x,y,z)**2
  enddo; enddo; enddo
!   ion kinetic energy
  ekin=0. ; ekinx=0. ; ekiny=0. ; ekinz=0.
  do i=1,np_i
    ekinx=ekinx+(rvi(4,i)**2)/(2*n_avg)
    ekiny=ekiny+(rvi(5,i)**2)/(2*n_avg)
    ekinz=ekinz+(rvi(6,i)**2)/(2*n_avg)
    ekin=ekin+(rvi(4,i)**2+rvi(5,i)**2+rvi(6,i)**2)/(2*n_avg)
  enddo
!   total energy density
  eelc=totalsum(eelc)*(dx*dy*dz)/(lx*ly*lz)
  eprs=totalsum(eprs)*(dx*dy*dz)/(lx*ly*lz)
  emag=totalsum(emag)*(dx*dy*dz)/(lx*ly*lz)
  emagx=totalsum(emagx)*(dx*dy*dz)/(lx*ly*lz)
  emagy=totalsum(emagy)*(dx*dy*dz)/(lx*ly*lz)
  emagz=totalsum(emagz)*(dx*dy*dz)/(lx*ly*lz)
  etaloss=totalsum(etaloss)*(dx*dy*dz)/(lx*ly*lz)
  

! Added by Tulasi to calculate the ion fluid current here
    eifl=0. ; eiflz=0. ; eiflx=0. ; eifly=0.
  do z=1,nz; do y=1,ny; do x=1,nx
    eiflx=eiflx+(jx(x,y,z)**2)/(2.*rho(x,y,z))
    eifly=eifly+(jy(x,y,z)**2)/(2.*rho(x,y,z))
    eiflz=eiflz+(jz(x,y,z)**2)/(2.*rho(x,y,z))
    eifl=eifl+(jx(x,y,z)**2+jy(x,y,z)**2+jz(x,y,z)**2)/(2.*rho(x,y,z))
  enddo; enddo; enddo

  do z=1,nz; do y=1,ny; do x=1,nx
    veloix(x,y,z)=jx(x,y,z)/rho(x,y,z)
    veloiy(x,y,z)=jy(x,y,z)/rho(x,y,z)
    veloiz(x,y,z)=jz(x,y,z)/rho(x,y,z)
  enddo; enddo; enddo

! Calculating the voritcities and enstrophies
b_enstrophy_x=0. ; b_enstrophy_y=0. ; b_enstrophy_z=0. 
v_enstrophy_x=0. ; v_enstrophy_y=0. ; v_enstrophy_z=0.
  do z=1,nz; do y=1,ny; do x=1,nx

    bomegax(x,y,z)=DY(b1z)-DZ(b1y)
    bomegay(x,y,z)=DZ(b1x)-DX(b1z)
    bomegaz(x,y,z)=DX(b1y)-DY(b1x)

    vomegax(x,y,z)=DY(veloiz)-DZ(veloiy)
    vomegay(x,y,z)=DZ(veloix)-DX(veloiz)
    vomegaz(x,y,z)=DX(veloiy)-DY(veloix)

    b_enstrophy_x=b_enstrophy_x+bomegax(x,y,z)**2
    b_enstrophy_y=b_enstrophy_y+bomegay(x,y,z)**2
    b_enstrophy_z=b_enstrophy_z+bomegaz(x,y,z)**2

    v_enstrophy_x=v_enstrophy_x+vomegax(x,y,z)**2
    v_enstrophy_y=v_enstrophy_y+vomegay(x,y,z)**2
    v_enstrophy_z=v_enstrophy_z+vomegaz(x,y,z)**2

  enddo; enddo; enddo



    eem=totalsum(eem)*(dx*dy*dz)/(lx*ly*lz) 

    ekin=totalsum(ekin)*(dx*dy*dz)/(lx*ly*lz) 
    ekinx=totalsum(ekinx)*(dx*dy*dz)/(lx*ly*lz)
    ekiny=totalsum(ekiny)*(dx*dy*dz)/(lx*ly*lz)
    ekinz=totalsum(ekinz)*(dx*dy*dz)/(lx*ly*lz)

    eifl=totalsum(eifl)*(dx*dy*dz)/(lx*ly*lz)
    eiflx=totalsum(eiflx)*(dx*dy*dz)/(lx*ly*lz)
    eifly=totalsum(eifly)*(dx*dy*dz)/(lx*ly*lz)
    eiflz=totalsum(eiflz)*(dx*dy*dz)/(lx*ly*lz)

    b_enstrophy_p=totalsum(b_enstrophy_x+b_enstrophy_y)*(dx*dy*dz)/(lx*ly*lz)
    b_enstrophy_x=totalsum(b_enstrophy_x)*(dx*dy*dz)/(lx*ly*lz)
    b_enstrophy_y=totalsum(b_enstrophy_y)*(dx*dy*dz)/(lx*ly*lz)
    b_enstrophy_z=totalsum(b_enstrophy_z)*(dx*dy*dz)/(lx*ly*lz)

    v_enstrophy_p=totalsum(v_enstrophy_x+v_enstrophy_y)*(dx*dy*dz)/(lx*ly*lz)
    v_enstrophy_x=totalsum(v_enstrophy_x)*(dx*dy*dz)/(lx*ly*lz)
    v_enstrophy_y=totalsum(v_enstrophy_y)*(dx*dy*dz)/(lx*ly*lz)
    v_enstrophy_z=totalsum(v_enstrophy_z)*(dx*dy*dz)/(lx*ly*lz)

  eges=eem+ekin+elviscen
  eprsix=ekinx-eiflx
  eprsiy=ekiny-eifly
  eprsiz=ekinz-eiflz
  eprsi =ekin -eifl
#endif

#ifdef twofluid
! two-fluid version
!   magnetic field energy, electron kinetic energy, electron thermal energy
!   ion kinetic energy, ion thermal energy
!   1/2 B**2/2 + 1/2 d_e**2 j**2 + p_e/(gam_e-1) + 1/2 ji**2/n +p_i/(gam_i-1)
  eges=0.
  call bound_b(b1x,b1y,b1z)
  do z=1,nz; do y=1,ny; do x=1,nx
    jtx(x,y,z)=DY(b1z)-DZ(b1y)
    jty(x,y,z)=DZ(b1x)-DX(b1z)
    jtz(x,y,z)=DX(b1y)-DY(b1x)
    eges=eges+0.5*(b1x(x,y,z)**2+b1y(x,y,z)**2+b1z(x,y,z)**2)+ &
              0.5*d_e2*(jtx(x,y,z)**2+jty(x,y,z)**2+jtz(x,y,z)**2)+ &
              pe1(x,y,z)/(gam_e-1.)+ pi1(x,y,z)/(gam_i-1.)+ &
              0.5*(j1x(x,y,z)**2+j1y(x,y,z)**2+j1z(x,y,z)**2)/n1(x,y,z)
  enddo; enddo; enddo
  eges=totalsum(eges)*(dx*dy*dz)/(lx*ly*lz) 
#endif

#ifdef hybrid
  if (myproc == 0) then
    write(9,'(30E14.5)') t, eges, emag, emagx, emagy, emagz, eelc, ekin, ekinx, ekiny, ekinz,&
               eifl, eiflx, eifly, eiflz, eprs, eprsi, eprsix, eprsiy, eprsiz, &
               v_enstrophy_p, v_enstrophy_x, v_enstrophy_y, v_enstrophy_z, &
               b_enstrophy_p, b_enstrophy_x, b_enstrophy_y, b_enstrophy_z, elviscen, etaloss
  endif
#endif
  call stop_timer(t_energy,n_energy)
End Subroutine energy

!------------------------------------------------------------------------------
!                      load balance statistics
!------------------------------------------------------------------------------

#ifdef fullparticle
Subroutine load_balance()
  Use pe_env
  Use partfield
  Implicit None
  Integer :: minnp_i,maxnp_i,avgnp_i,minnp_e,maxnp_e,avgnp_e,pe
  Integer :: mpi_err

  call MPI_Allreduce(np_i,minnp_i,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(np_i,maxnp_i,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(np_i,avgnp_i,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(np_e,minnp_e,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(np_e,maxnp_e,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(np_e,avgnp_e,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  avgnp_i=avgnp_i/n_pes; avgnp_e=avgnp_e/n_pes

  if (myproc == 0) then
    write(6,*) '--- load balance statistics ---'
    write(6,*) 'ions per PE (avg/min/max):',avgnp_i,'/',minnp_i,'/',maxnp_i
    write(6,*) (maxnp_i*100.)/maxparticles,'% of ion buffer used'
    write(6,*) 'electrons per PE (avg/min/max):',avgnp_e,'/',minnp_e,'/',maxnp_e
    write(6,*) (maxnp_e*100.)/maxparticles,'% of electron buffer used'
    call flush(6)
  endif
End Subroutine load_balance
#endif

#ifdef hybrid
Subroutine load_balance()
  Use pe_env
  Use partfield
  Implicit None
  Integer :: minnp_i,maxnp_i,avgnp_i,pe
  Integer :: mpi_err

  call MPI_Allreduce(np_i,minnp_i,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(np_i,maxnp_i,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(np_i,avgnp_i,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  avgnp_i=avgnp_i/n_pes

  if (myproc == 0) then
    write(6,*) '--- load balance statistics ---'
    write(6,*) 'ions per PE (avg/min/max):',avgnp_i,'/',minnp_i,'/',maxnp_i
    write(6,*) (maxnp_i*100.)/maxparticles,'% of ion buffer used'
    call flush(6)
  endif
End Subroutine load_balance
#endif

!------------------------------------------------------------------------------
!                             subtract average
!------------------------------------------------------------------------------

Subroutine subtract_average(data)

  Use pe_env
  Implicit None
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: data
  Real(kind=drk) avg
  Integer :: x,y,z

  avg=0.
  do z=1,nz; do y=1,ny; do x=1,nx
    avg=avg+data(x,y,z)
  enddo; enddo; enddo
  avg=totalsum(avg)/(nx*ny*nz*pex*pey*pez)
  do z=1,nz; do y=1,ny; do x=1,nx
    data(x,y,z)=data(x,y,z)-avg
  enddo; enddo; enddo
End Subroutine subtract_average
