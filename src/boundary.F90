!******************************************************************************
!                3D particle code: non periodic boundary conditions
!                         Andreas Zeiler, 1999, 2000
!
!                         LATEST CHANGE: May 10, 2001
!
!******************************************************************************

#include "param"

! assign unique number to each set of boundary conditions
#define periodic 0
#define gem_reconnection 1
#define gem_3d 2

!------------------------------------------------------------------------------
!        boundary conditions for electric and magnetic field,
!        potential phi, current density, charge density, electron pressure
!        ion pressure, density
!------------------------------------------------------------------------------

Subroutine bound_e(ex,ey,ez)
! boundary conditions for the electric field
  Use pe_env
  Implicit None
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: ex,ey,ez

  call boundary(0,0,ex); call boundary(0,0,ey); call boundary(0,0,ez)
#if (boundary_condition == periodic)
! periodic boundary conditions (nothing to be done)
#elif (boundary_condition == gem_reconnection)
! gem reconnection setup of harris sheet (Mike Shay)
  call von_neumann(0,0,ex)
  call dirichlet(0,0,ey)
  call dirichlet(0,0,ez)
#elif (boundary_condition == gem_3d)
! gem reconnection setup of harris sheet, 3D case with conducting wall
  call dirichlet(0,0,ex)
  call von_neumann(0,0,ey)
  call dirichlet(0,0,ez)
#else
  we never should get here
#endif
End Subroutine bound_e

!------------------------------------------------------------------------------

Subroutine bound_b(bx,by,bz)
! boundary conditions for the magnetic field
  Use pe_env
  Implicit None
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: bx,by,bz

  call boundary(0,0,bx); call boundary(0,0,by); call boundary(0,0,bz)
#if (boundary_condition == periodic)
! periodic boundary conditions (nothing to be done)
#elif (boundary_condition == gem_reconnection)
! gem reconnection setup of harris sheet (Mike Shay)
  call von_neumann(0,0,bx)
  call dirichlet(0,0,by)
  call dirichlet(0,0,bz)
#elif (boundary_condition == gem_3d)
! gem reconnection setup of harris sheet, 3D case with conducting wall
  call von_neumann(0,0,bx)
  call dirichlet(0,0,by)
  call von_neumann(0,0,bz)
#else
  we never should get here
#endif
End Subroutine bound_b

!------------------------------------------------------------------------------

Subroutine bound_phi(levelxy,levelz,phi)
! boundary conditions for potential phi
  Use pe_env
  Implicit None
  Integer :: levelxy,levelz
  Real(kind=drk), Dimension(0:nx/2**levelxy+1,0:ny/2**levelxy+1,0:nz/2**levelz+1) :: phi

  call boundary(levelxy,levelz,phi)
#if (boundary_condition == periodic)
! periodic boundary conditions (nothing to be done)
#elif (boundary_condition == gem_reconnection)
! gem reconnection setup of harris sheet (Mike Shay)
  call von_neumann(levelxy,levelz,phi)
#elif (boundary_condition == gem_3d)
! gem reconnection setup of harris sheet, 3D case with conducting wall
  call dirichlet(levelxy,levelz,phi)
#else
  we never should get here
#endif
End Subroutine bound_phi

!------------------------------------------------------------------------------

Subroutine bound_j(jx,jy,jz)
! boundary conditions for the current density
  Use pe_env
  Implicit None
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: jx,jy,jz

  call boundary(0,0,jx); call boundary(0,0,jy); call boundary(0,0,jz)
#if (boundary_condition == periodic)
! periodic boundary conditions (nothing to be done)
#elif (boundary_condition == gem_reconnection)
! gem reconnection setup of harris sheet (Mike Shay)
  call von_neumann(0,0,jx)
  call dirichlet(0,0,jy)
  call dirichlet(0,0,jz)
#elif (boundary_condition == gem_3d)
! gem reconnection setup of harris sheet, 3D case with conducting wall
  call dirichlet(0,0,jx)
  call von_neumann(0,0,jy)
  call dirichlet(0,0,jz)
#else
  we never should get here
#endif
End Subroutine bound_j

!------------------------------------------------------------------------------

Subroutine bound_rho(rho)
! boundary conditions for the charge density
  Use pe_env
  Implicit None
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: rho

  call boundary(0,0,rho)
#if (boundary_condition == periodic)
! periodic boundary conditions (nothing to be done)
#elif (boundary_condition == gem_reconnection)
! gem reconnection setup of harris sheet (Mike Shay)
  call von_neumann(0,0,rho)
#elif (boundary_condition == gem_3d)
! gem reconnection setup of harris sheet, 3D case with conducting wall
  call dirichlet(0,0,rho)
#else
  we never should get here
#endif
End Subroutine bound_rho

!------------------------------------------------------------------------------

Subroutine bound_pe(pe)
! boundary conditions for the electron pressure
! (hybrid and two-fluid version)
  Use pe_env
  Implicit None
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: pe

  call boundary(0,0,pe)
#if (boundary_condition == periodic)
! periodic boundary conditions (nothing to be done)
#else
  if (myproc==0) then
    write(6,*) 'only periodic boundary conditions defined for electron pressure'
    call exitallpes()
  endif
#endif
End Subroutine bound_pe

!------------------------------------------------------------------------------

Subroutine bound_pi(pi)
! boundary conditions for the ion pressure
! (two-fluid version)
  Use pe_env
  Implicit None
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: pi

  call boundary(0,0,pi)
#if (boundary_condition == periodic)
! periodic boundary conditions (nothing to be done)
#else
  if (myproc==0) then
    write(6,*) 'only periodic boundary conditions defined for ion pressure'
    call exitallpes()
  endif
#endif
End Subroutine bound_pi

!------------------------------------------------------------------------------

Subroutine bound_n(n)
! boundary conditions for the fluid density
! (two-fluid version)
  Use pe_env
  Implicit None
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: n

  call boundary(0,0,n)
#if (boundary_condition == periodic)
! periodic boundary conditions (nothing to be done)
#else
  if (myproc==0) then
    write(6,*) 'only periodic boundary conditions defined for fluid density'
    call exitallpes()
  endif
#endif
End Subroutine bound_n

!------------------------------------------------------------------------------
!                     virtual charges and virtual currents
!------------------------------------------------------------------------------

Subroutine virtual_charge(mc)
! return virtual charge (called by calc_rho)
  Use pe_env
  Implicit None
  Integer :: mc
 
#if (boundary_condition == periodic)
! periodic boundary conditions: no virtual particles
  mc=0
#elif (boundary_condition == gem_reconnection)
! gem reconnection setup of harris sheet (Mike Shay)
  mc=1                             ! virtual charge with same sign
#elif (boundary_condition == gem_3d)
! gem reconnection setup of harris sheet, 3D case with conducting wall
  mc=-1                            ! virtual charge with opposite sign
#else
  we never should get here
#endif
End Subroutine virtual_charge

!------------------------------------------------------------------------------

Subroutine virtual_current(mcx,mcy,mcz)
! return virtual current (called by calc_j)
  Use pe_env
  Implicit None
  Integer :: mcx,mcy,mcz
 
#if (boundary_condition == periodic)
! periodic boundary conditions: no virtual particles
  mcx=0; mcy=0; mcz=0
#elif (boundary_condition == gem_reconnection)
! gem reconnection setup of harris sheet (Mike Shay)
  mcx = 1                         ! virtual particle: same momentum 
  mcy =-1                         ! in x, opposite momentum in y and z   
  mcz =-1
#elif (boundary_condition == gem_3d)
! gem reconnection setup of harris sheet, 3D case with conducting wall
  mcx = 1                         ! virtual particle: same momentum 
  mcy =-1                         ! in x and z, opposite momentum in y   
  mcz = 1
#else
  we never should get here
#endif
End Subroutine virtual_current

!------------------------------------------------------------------------------

Function non_periodic()
! switch for transfer of boundary layers in bound_rho_j and for routine stepx
  Logical :: non_periodic

#if (boundary_condition == periodic)
! periodic boundary conditions
  non_periodic=.false.
#elif (boundary_condition == gem_reconnection)
! gem reconnection setup of harris sheet (Mike Shay)
  non_periodic=.true.
#elif (boundary_condition == gem_3d)
! gem reconnection setup of harris sheet, 3D case with conducting wall
  non_periodic=.true.
#else
  we never should get here
#endif
  return
End Function non_periodic

!------------------------------------------------------------------------------
!                     particles at y boundary
!------------------------------------------------------------------------------

Subroutine bound_particle(ix,iy,iz)
! specifies which components of the particle momentum need to be reverted
! when particle hits boundary in y
  Use pe_env
  Use partfield
  Implicit None
  Integer :: ix,iy,iz

#if (boundary_condition == periodic)
! periodic boundary conditions (nothing to be done)
#elif (boundary_condition == gem_reconnection)
! gem reconnection setup of harris sheet (Mike Shay)
  ix = 1                             ! leave momentum in x unchanged
  iy =-1                             ! invert momentum in y and z
  iz =-1
#elif (boundary_condition == gem_3d)
! gem reconnection setup of harris sheet, 3D case with conducting wall
  ix = 1                             ! leave momentum in x and z unchanged
  iy =-1                             ! invert momentum in y
  iz = 1
#else
  we never should get here
#endif
End Subroutine bound_particle

!------------------------------------------------------------------------------
!          3D-boundary condition (communication between subdomains)
!------------------------------------------------------------------------------

Subroutine boundary(levelxy,levelz,f)
! exchanges the boundary layers between the different volumes (PEs) in the
! domain decomposition; periodic boundary conditions in all directions

  Use pe_env ! Contains also mpif.h
  Implicit None

  Integer levelxy,levelz,nnx,nny,nnz,icnt
  Real(kind=drk), Dimension(0:nx/2**levelxy+1,0:ny/2**levelxy+1, &
                            0:nz/2**levelz+1) :: f
  Real(kind=drk), Dimension(1:ny/2**levelxy  ,1:nz/2**levelz  )  :: tmpxs, tmpxr
  Real(kind=drk), Dimension(0:nx/2**levelxy+1,1:nz/2**levelz  )  :: tmpys, tmpyr
  Real(kind=drk), Dimension(0:nx/2**levelxy+1,0:ny/2**levelxy+1) :: tmpzs, tmpzr
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status
  Integer :: mpi_err

  nnx=nx/2**levelxy; nny=ny/2**levelxy; nnz=nz/2**levelz

! Exchange in X-Direction

  icnt = nny*nnz

  tmpxs(1:nny,1:nnz) = f(1,1:nny,1:nnz)
  call MPI_Sendrecv(tmpxs,icnt,mpi_drk,nb_pe(-1,0,0),111, &
                    tmpxr,icnt,mpi_drk,nb_pe( 1,0,0),111, &
                    MPI_COMM_WORLD,mpi_status,mpi_err)
  f(nnx+1,1:nny,1:nnz) = tmpxr(1:nny,1:nnz)

  tmpxs(1:nny,1:nnz) = f(nnx,1:nny,1:nnz)
  call MPI_Sendrecv(tmpxs,icnt,mpi_drk,nb_pe( 1,0,0),222, &
                    tmpxr,icnt,mpi_drk,nb_pe(-1,0,0),222, &
                    MPI_COMM_WORLD,mpi_status,mpi_err)
  f(0,1:nny,1:nnz) = tmpxr(1:nny,1:nnz)

! Exchange in Y-Direction

  icnt = (nnx+2)*nnz

  tmpys(0:nnx+1,1:nnz) = f(0:nnx+1,1,1:nnz)
  call MPI_Sendrecv(tmpys,icnt,mpi_drk,nb_pe(0,-1,0),333, &
                    tmpyr,icnt,mpi_drk,nb_pe(0, 1,0),333, &
                    MPI_COMM_WORLD,mpi_status,mpi_err)
  f(0:nnx+1,nny+1,1:nnz) = tmpyr(0:nnx+1,1:nnz)

  tmpys(0:nnx+1,1:nnz) = f(0:nnx+1,nny,1:nnz)
  call MPI_Sendrecv(tmpys,icnt,mpi_drk,nb_pe(0, 1,0),444, &
                    tmpyr,icnt,mpi_drk,nb_pe(0,-1,0),444, &
                    MPI_COMM_WORLD,mpi_status,mpi_err)
  f(0:nnx+1,0,1:nnz) = tmpyr(0:nnx+1,1:nnz)

! Exchange in Z-Direction
  if (exchange_z_boundary) then
    icnt = (nnx+2)*(nny+2)

    tmpzs(0:nnx+1,0:nny+1) = f(0:nnx+1,0:nny+1,1)
    call MPI_Sendrecv(tmpzs,icnt,mpi_drk,nb_pe(0,0,-1),555, &
                      tmpzr,icnt,mpi_drk,nb_pe(0,0, 1),555, &
                      MPI_COMM_WORLD,mpi_status,mpi_err)
    f(0:nnx+1,0:nny+1,nnz+1) = tmpzr(0:nnx+1,0:nny+1)

    tmpzs(0:nnx+1,0:nny+1) = f(0:nnx+1,0:nny+1,nnz)
    call MPI_Sendrecv(tmpzs,icnt,mpi_drk,nb_pe(0,0, 1),666, &
                      tmpzr,icnt,mpi_drk,nb_pe(0,0,-1),666, &
                      MPI_COMM_WORLD,mpi_status,mpi_err)
    f(0:nnx+1,0:nny+1,0) = tmpzr(0:nnx+1,0:nny+1)
  endif

  return
End Subroutine boundary

!------------------------------------------------------------------------------
!               Dirichlet, von Neumann boundary conditions in y
!------------------------------------------------------------------------------

Subroutine dirichlet(levelxy,levelz,f)
! apply Dirichlet boundary conditions (f=0)
  Use pe_env
  Implicit none

  Integer :: levelxy,levelz,nnx,nny,nnz,x,z
  Real(kind=drk), Dimension(0:nx/2**levelxy+1,0:ny/2**levelxy+1,0:nz/2**levelz+1) :: f

  nnx=nx/2**levelxy; nny=ny/2**levelxy; nnz=nz/2**levelz
  if (my_pey == 0) then
    do z=0,nnz+1; do x=0,nnx+1
      f(x,0,z)=-f(x,1,z)
    enddo; enddo
  endif
  if (my_pey == pey-1) then
    do z=0,nnz+1; do x=0,nnx+1
      f(x,nny+1,z)=-f(x,nny,z)
    enddo; enddo
  endif
End Subroutine dirichlet

!------------------------------------------------------------------------------

Subroutine von_neumann(levelxy,levelz,f)
! apply von Neumann boundary conditions (df/dy=0)
  Use pe_env
  Implicit none

  Integer :: levelxy,levelz,nnx,nny,nnz,x,z
  Real(kind=drk), Dimension(0:nx/2**levelxy+1,0:ny/2**levelxy+1,0:nz/2**levelz+1) :: f

  nnx=nx/2**levelxy; nny=ny/2**levelxy; nnz=nz/2**levelz
  if (my_pey == 0) then
    do z=0,nnz+1; do x=0,nnx+1
      f(x,0,z)=f(x,1,z)
    enddo; enddo
  endif
  if (my_pey == pey-1) then
    do z=0,nnz+1; do x=0,nnx+1
      f(x,nny+1,z)=f(x,nny,z)
    enddo; enddo
  endif
End Subroutine von_neumann

!------------------------------------------------------------------------------
!                smooth electric field, current and charge 
!------------------------------------------------------------------------------

Subroutine smooth_e(ex,ey,ez,esx,esy,esz)
! smooth electric field
  Use pe_env
  Implicit None
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: ex,ey,ez,esx,esy,esz

  call bound_e(ex,ey,ez)
  call smooth(ex,esx); call smooth(ey,esy); call smooth(ez,esz)
End Subroutine smooth_e

!------------------------------------------------------------------------------

Subroutine smooth_b(bx,by,bz,bsx,bsy,bsz)
! smooth magnetic field
  Use pe_env
  Implicit None
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: bx,by,bz,bsx,bsy,bsz

  call bound_b(bx,by,bz)
  call smooth(bx,bsx); call smooth(by,bsy); call smooth(bz,bsz)
End Subroutine smooth_b

!------------------------------------------------------------------------------

Subroutine smooth_j(jx,jy,jz,jsx,jsy,jsz)
! smooth current density
  Use pe_env
  Use timer
  Implicit None
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: jx,jy,jz,jsx,jsy,jsz

  call start_timer()
  call bound_j(jx,jy,jz)
  call smooth(jx,jsx); call smooth(jy,jsy); call smooth(jz,jsz)
  call stop_timer(t_smooth_j,n_smooth_j)
End Subroutine smooth_j

!------------------------------------------------------------------------------

Subroutine smooth_rho(rho,rhos)
! smooth charge density
  Use pe_env
  Use timer
  Implicit None
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: rho,rhos

  call start_timer()
  call bound_rho(rho)
  call smooth(rho,rhos)
  call stop_timer(t_smooth_rho,n_smooth_rho)
End Subroutine smooth_rho

!------------------------------------------------------------------------------

Subroutine smooth(input,output)
! smooth scalar field
  Use pe_env
  Implicit None

! To be consistent, this smoothing function should be changed when the
! fourth_order flag is set.  With the usual derivative (second-order),
! the smoothing and differentiation operators are commutative.  I do not
! believe that is generally true.

  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: input, output, tmp
  Integer :: x,y,z

! boundary layers of input array must be set appropriately by calling routine
  do z=1,nz; do y=1,ny; do x=1,nx
    tmp(x,y,z)= 1/64.* ( input(x+1,y+1,z+1) + input(x-1,y+1,z+1) + &
                         input(x+1,y-1,z+1) + input(x-1,y-1,z+1) + &
                         input(x+1,y+1,z-1) + input(x-1,y+1,z-1) + &
                         input(x+1,y-1,z-1) + input(x-1,y-1,z-1)   ) + &
                1/32.* ( input(x  ,y+1,z+1) + input(x  ,y-1,z+1) + &
                         input(x  ,y+1,z-1) + input(x  ,y-1,z-1) + &
                         input(x+1,y  ,z+1) + input(x-1,y  ,z+1) + &
                         input(x+1,y  ,z-1) + input(x-1,y  ,z-1) + &
                         input(x+1,y+1,z  ) + input(x-1,y+1,z  ) + &
                         input(x+1,y-1,z  ) + input(x-1,y-1,z  ) ) + &
                1/16.* ( input(x+1,y  ,z  ) + input(x-1,y  ,z  ) + &
                         input(x  ,y+1,z  ) + input(x  ,y-1,z  ) + &
                         input(x  ,y  ,z+1) + input(x  ,y  ,z-1)   ) + &
                1/8. *   input(x  ,y  ,z  )
  enddo; enddo; enddo
  do z=1,nz; do y=1,ny; do x=1,nx
    output(x,y,z)=tmp(x,y,z)
  enddo; enddo; enddo
End Subroutine smooth











