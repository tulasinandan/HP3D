!******************************************************************************
!                   3D particle code: Multigrid-Solver
!                           Andreas Zeiler, 1999
!			    Marc Swisdak, 2004
!
!                       LATEST CHANGE: June 23, 2004
!
!******************************************************************************

#include "param"
#ifndef miniter
#define miniter 0
#endif
#ifndef fixiter
#define fixiter 6
#endif
#ifndef c_2
#define c_2 1.
#endif
#ifndef d_e2
#define d_e2 1.
#endif

#if (nz == 1 && pez == 1)
#define iz 0
#else
#define iz 1
#endif

#ifndef GS_LEX
#ifndef GS_RB
#ifndef JACOBI
#define JACOBI
#endif
#endif
#endif

Module mgctrl
  Implicit none
  Real :: c1, c2,c3, normalize
  Logical :: use_sor
End Module mgctrl

!------------------------------------------------------------------------------
!                      call as poisson solver
!------------------------------------------------------------------------------

Subroutine poisson(f,b,converg)
! solve equation delt**2 f=b by multigrid iteration
! iteration stops when residual < converg and iters >= miniter

  Use pe_env
  Use mgctrl
  Implicit None
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: f,b
  Real(kind=drk) :: converg

! set constants for Poisson equation
  c1=0; c2=1; c3=-1
! set normalization constant for residual
#ifdef norm_res
  normalize=nx*pex*ny*pey*nz*pez*c_2
#else
  normalize=1.
#endif
! switch-on successive overrelaxation
#ifdef coarsegrid_sor
  use_sor=.true.
#else
  use_sor=.false.
#endif
  call multigrid(f,b,converg,.false.)
End Subroutine poisson

!------------------------------------------------------------------------------
!                      call as helmholtz solver
!------------------------------------------------------------------------------

Subroutine helmholtz(f,b,converg,fixed_cycles)
! solve equation (1-de2*delt**2) f=b by multigrid iteration
! iteration stops when residual < converg and iters >= miniter

  Use pe_env
  Use mgctrl
  Implicit None
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: f,b
  Real(kind=drk) :: converg
  Logical :: fixed_cycles

! set constants for Helmholtz equation
  c1=1; c2=d_e2; c3=1
! set normalization constant for residual
#ifdef norm_res
  normalize=nx*pex*ny*pey*nz*pez
#else
  normalize=1.
#endif
! no successive overrelaxation in helmholtz case
  use_sor=.false.
  call multigrid(f,b,converg,fixed_cycles)
End Subroutine helmholtz

!------------------------------------------------------------------------------
!                      3D-multigrid / main entry point
!------------------------------------------------------------------------------

Subroutine multigrid(f,b,converg,fixed_cycles)
! solve equation (c1-c2*delt**2) f=c3*b
! iteration stops when residual < converg and iters >= miniter

  Use pe_env
  Implicit None
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: f,b
  Real(kind=drk) :: converg, res3d, res
  Integer :: iter,pe,x,y,z
  Logical :: fixed_cycles 

! average charge must be zero, otherwise multigrid fails to converge
  if (fixed_cycles) then
    do iter=1,fixiter
      if (nz .ge. 2) then
        call mg3d_0(f,b)
      else
        if (pez .eq. 1) exchange_z_boundary=.false.
        call mg3dxy_0(f,b)
        exchange_z_boundary=.true.
      endif
    enddo
#ifdef mgverbose
    if (myproc == 0) write(6,*) 'multigrid: ',fixiter,' iterations (fixed)'
#endif
  else
    res=res3d(f,b)
#ifdef mgverbose2
      if (myproc==0) write(6,*) 'multigrid: ',res
#endif
    iter=0
    do while (res .gt. converg .or. iter .lt. miniter)
      if (iter .gt. maxiter) then
        if (myproc == 0) write(6,*) '***** multigrid: failed to converge *****'
        call exitallpes()
      endif
      if (nz .ge. 2) then
        call mg3d_0(f,b)
      else
        if (pez .eq. 1) exchange_z_boundary=.false.
        call mg3dxy_0(f,b)
        exchange_z_boundary=.true.
      endif
      res=res3d(f,b)
#ifdef mgverbose2
      if (myproc==0) write(6,*) 'multigrid: ',res
#endif
      iter=iter+1
    enddo
#ifdef mgverbose
    if (myproc == 0) write(6,*) 'multigrid: ',iter,' iterations'
#endif
  endif
  return
End Subroutine multigrid

!------------------------------------------------------------------------------
!                          3D-multigrid hierarchy
!------------------------------------------------------------------------------

! preprocessor generates subroutine (level xx) for each grid level

! preprocessor command: expand operator into finite differences
define(DELTA,`&
  (($1(x+1,y,z)-2*$1(x,y,z)+$1(x-1,y,z))/dx**2+&
   ($1(x,y+1,z)-2*$1(x,y,z)+$1(x,y-1,z))/dy**2+&
   ($1(x,y,z+iz)-2*$1(x,y,z)+$1(x,y,z-iz))/dz**2)&
  ')
define(DELTA2,`&
  (($1(x+1,y,1)-2*$1(x,y,1)+$1(x-1,y,1))/dx**2+&
   ($1(x,y+1,1)-2*$1(x,y,1)+$1(x,y-1,1))/dy**2+&
   ($1(x,y,1+iz)-2*$1(x,y,1)+$1(x,y,1-iz))/dz**2)&
  ')

! preprocessor command: generate subgrid at level xx
define(CREATE_MG,`
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                          3D-multigrid / level $1
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Subroutine mg3d_$1(f,b)
! call if nz >=2 at current grid level
  Use pe_env
  Use mgctrl
  Implicit None

  Real(kind=drk), Dimension(0:nx/2**$1+1,0:ny/2**$1+1,0:nz/2**$1+1) :: f,b
  Real(kind=drk), Save, Dimension(0:nx/2**$1+1,0:ny/2**$1+1,0:nz/2**$1+1) :: g
  Real(kind=drk), Save, &
     Dimension(0:nx/2**($1+1)+1,0:ny/2**($1+1)+1,0:nz/2**($1+1)+1) :: fc,bc
  Integer :: x,y,z,i,j,nu4
  Real(kind=drk), Parameter :: dx=lx/(nx*pex)*2**$1, dy=ly/(ny*pey)*2**$1, &
      dz=lz/(nz*pez)*2**$1
  Real(kind=drk) :: dtau

! calculate time step dtau for smoother
! cpp commands added to deal with possibility that dz is arbitrary in 2D
! simulation
#if (nz == 1 && pez == 1)
   dtau=1./(c2*(4./dx**2+4./dy**2)+c1)
#else
   dtau=1./(c2*(4./dx**2+4./dy**2+4./dz**2)+c1)
#endif


! pre-smoothing (nu1 iterationsi, choose JACOBI, GS_LEX, GS_RB)
#ifdef JACOBI
  do i=1,nu1
    call bound_phi($1,$1,f)
    do z=1,nz/2**$1; do y=1,ny/2**$1; do x=1,nx/2**$1
      g(x,y,z)=c2*DELTA(f)-c1*f(x,y,z)+c3*b(x,y,z)
    enddo; enddo; enddo
    do z=1,nz/2**$1; do y=1,ny/2**$1; do x=1,nx/2**$1
      f(x,y,z)=f(x,y,z)+dtau*g(x,y,z)
    enddo; enddo; enddo
  enddo
#endif
#ifdef GS_LEX
  do i=1,nu1
    call bound_phi($1,$1,f)
    do z=1,nz/2**$1; do y=1,ny/2**$1; do x=1,nx/2**$1
      f(x,y,z)=f(x,y,z)+dtau*(c2*DELTA(f)-c1*f(x,y,z)+c3*b(x,y,z))
    enddo; enddo; enddo
  enddo
#endif
#ifdef GS_RB
  do i=1,2*nu1
    call bound_phi($1,$1,f)
    do z=1,nz/2**$1; do y=1,ny/2**$1; do x=1,nx/2**$1
      if (mod(x+my_pex*nx/2**$1+y+my_pey*ny/2**$1+ &
              z+my_pez*nz/2**$1,2)==mod(i,2)) &
      f(x,y,z)=f(x,y,z)+dtau*(c2*DELTA(f)-c1*f(x,y,z)+c3*b(x,y,z))
    enddo; enddo; enddo
  enddo
#endif
! multigrid recursion
  if (mod(nx/2**$1,2)+mod(ny/2**$1,2) == 0) then
! restriction of defect (27 or 8 point stencil), guess for coarse grid solution
    call bound_phi($1,$1,f)
    do z=1,nz/2**$1; do y=1,ny/2**$1; do x=1,nx/2**$1
      g(x,y,z)=-c3*(c2*DELTA(f)-c1*f(x,y,z)+c3*b(x,y,z))
    enddo; enddo; enddo
#ifdef alt_restrict_prolong
    call boundary($1,$1,g)
    do z=1,nz/2**($1+1); do y=1,ny/2**($1+1); do x=1,nx/2**($1+1)
      bc(x,y,z)=(g(2*x  ,2*y  ,2*z  )   +&
                 g(2*x-1,2*y  ,2*z  )/2 + g(2*x+1,2*y  ,2*z  )/2 +&
                 g(2*x  ,2*y-1,2*z  )/2 + g(2*x  ,2*y+1,2*z  )/2 +&
                 g(2*x  ,2*y  ,2*z-1)/2 + g(2*x  ,2*y  ,2*z+1)/2 +&
                 g(2*x-1,2*y-1,2*z  )/4 + g(2*x-1,2*y+1,2*z  )/4 +&
                 g(2*x+1,2*y-1,2*z  )/4 + g(2*x+1,2*y+1,2*z  )/4 +&
                 g(2*x-1,2*y  ,2*z-1)/4 + g(2*x-1,2*y  ,2*z+1)/4 +&
                 g(2*x+1,2*y  ,2*z-1)/4 + g(2*x+1,2*y  ,2*z+1)/4 +&
                 g(2*x  ,2*y-1,2*z-1)/4 + g(2*x  ,2*y-1,2*z+1)/4 +&
                 g(2*x  ,2*y+1,2*z-1)/4 + g(2*x  ,2*y+1,2*z+1)/4 +&
                 g(2*x-1,2*y-1,2*z-1)/8 + g(2*x-1,2*y-1,2*z+1)/8 +&
                 g(2*x-1,2*y+1,2*z-1)/8 + g(2*x+1,2*y-1,2*z-1)/8 +&
                 g(2*x-1,2*y+1,2*z+1)/8 + g(2*x+1,2*y-1,2*z+1)/8 +&
                 g(2*x+1,2*y+1,2*z-1)/8 + g(2*x+1,2*y+1,2*z+1)/8   )/8
      fc(x,y,z)=0.
    enddo; enddo; enddo
#else
    do z=1,nz/2**($1+1); do y=1,ny/2**($1+1); do x=1,nx/2**($1+1)
      bc(x,y,z)=(g(2*x  ,2*y  ,2*z  ) + g(2*x-1,2*y  ,2*z  ) +&
                 g(2*x  ,2*y-1,2*z  ) + g(2*x-1,2*y-1,2*z  ) +&
                 g(2*x  ,2*y  ,2*z-1) + g(2*x-1,2*y  ,2*z-1) +&
                 g(2*x  ,2*y-1,2*z-1) + g(2*x-1,2*y-1,2*z-1))/8
      fc(x,y,z)=0.
    enddo; enddo; enddo
#endif
! multigrid recursion
ifelse($1,$2,,`
    if (nz/2**($1+1) .ge. 2) then
      call mg3d_`'eval($1+1)(fc,bc)
    else
      call mg3dxy_`'eval($1+1)(fc,bc)
    endif
')
#ifdef WCYCLE
ifelse($1,$2,,`
    if (nz/2**($1+1) .ge. 2) then
      call mg3d_`'eval($1+1)(fc,bc)
    else
      call mg3dxy_`'eval($1+1)(fc,bc)
    endif
')
#endif
! subtract prolongation (27 or 8 point stencil)
#ifdef alt_restrict_prolong
    call boundary($1+1,$1+1,fc)
    do z=1,nz/2**($1+1); do y=1,ny/2**($1+1); do x=1,nx/2**($1+1)
      f(2*x  ,2*y  ,2*z  ) = f(2*x  ,2*y  ,2*z  ) - fc(x,y,z)
      f(2*x-1,2*y  ,2*z  ) = f(2*x-1,2*y  ,2*z  ) &
         -(fc(x,y,z)+fc(x-1,y,z))/2
      f(2*x  ,2*y-1,2*z  ) = f(2*x  ,2*y-1,2*z  ) &
         -(fc(x,y,z)+fc(x,y-1,z))/2
      f(2*x  ,2*y  ,2*z-1) = f(2*x  ,2*y  ,2*z-1) &
         -(fc(x,y,z)+fc(x,y,z-1))/2
      f(2*x-1,2*y-1,2*z  ) = f(2*x-1,2*y-1,2*z  ) &
           -(fc(x,y,z)+fc(x,y-1,z)+fc(x-1,y,z)+fc(x-1,y-1,z))/4
      f(2*x-1,2*y  ,2*z-1) = f(2*x-1,2*y  ,2*z-1) &
           -(fc(x,y,z)+fc(x,y,z-1)+fc(x-1,y,z)+fc(x-1,y,z-1))/4
      f(2*x  ,2*y-1,2*z-1) = f(2*x  ,2*y-1,2*z-1) &
           -(fc(x,y,z)+fc(x,y,z-1)+fc(x,y-1,z)+fc(x,y-1,z-1))/4
      f(2*x-1,2*y-1,2*z-1) = f(2*x-1,2*y-1,2*z-1) &
           -(fc(x,y,z)+fc(x,y-1,z)+fc(x-1,y,z)+fc(x-1,y-1,z)+ &
             fc(x,y,z-1)+fc(x,y-1,z-1)+fc(x-1,y,z-1)+fc(x-1,y-1,z-1))/8
    enddo; enddo; enddo
#else
    do z=1,nz/2**($1+1); do y=1,ny/2**($1+1); do x=1,nx/2**($1+1)
      f(2*x  ,2*y  ,2*z  ) = f(2*x  ,2*y  ,2*z  ) - fc(x,y,z)
      f(2*x-1,2*y  ,2*z  ) = f(2*x-1,2*y  ,2*z  ) - fc(x,y,z)
      f(2*x  ,2*y-1,2*z  ) = f(2*x  ,2*y-1,2*z  ) - fc(x,y,z)
      f(2*x-1,2*y-1,2*z  ) = f(2*x-1,2*y-1,2*z  ) - fc(x,y,z)
      f(2*x  ,2*y  ,2*z-1) = f(2*x  ,2*y  ,2*z-1) - fc(x,y,z)
      f(2*x-1,2*y  ,2*z-1) = f(2*x-1,2*y  ,2*z-1) - fc(x,y,z)
      f(2*x  ,2*y-1,2*z-1) = f(2*x  ,2*y-1,2*z-1) - fc(x,y,z)
      f(2*x-1,2*y-1,2*z-1) = f(2*x-1,2*y-1,2*z-1) - fc(x,y,z)
    enddo; enddo; enddo
#endif
    nu4=nu2
  else
! additional smoothing iterations at coarsest grid
    nu4=nu2+nu3
    if (use_sor) then
! successive overrelaxation at coarsest grid
      call sor($1,$1,f,b)
      nu4=0
    endif
  endif
! post-smoothing (nu2 or nu2+nu3 iterations, choose JACOBI, GS_LEX, GS_RB)
#ifdef JACOBI
  do i=1,nu4
    call bound_phi($1,$1,f)
    do z=1,nz/2**$1; do y=1,ny/2**$1; do x=1,nx/2**$1
      g(x,y,z)=c2*DELTA(f)-c1*f(x,y,z)+c3*b(x,y,z)
    enddo; enddo; enddo
    do z=1,nz/2**$1; do y=1,ny/2**$1; do x=1,nx/2**$1
      f(x,y,z)=f(x,y,z)+dtau*g(x,y,z)
    enddo; enddo; enddo
  enddo
#endif
#ifdef GS_LEX
  do i=1,nu4
    call bound_phi($1,$1,f)
    do z=1,nz/2**$1; do y=1,ny/2**$1; do x=1,nx/2**$1
      f(x,y,z)=f(x,y,z)+dtau*(c2*DELTA(f)-c1*f(x,y,z)+c3*b(x,y,z))
    enddo; enddo; enddo
  enddo
#endif
#ifdef GS_RB
  do i=1,2*nu4
    call bound_phi($1,$1,f)
    do z=1,nz/2**$1; do y=1,ny/2**$1; do x=1,nx/2**$1
      if (mod(x+my_pex*nx/2**$1+y+my_pey*ny/2**$1+ &
              z+my_pez*nz/2**$1,2)==mod(i,2)) &
      f(x,y,z)=f(x,y,z)+dtau*(c2*DELTA(f)-c1*f(x,y,z)+c3*b(x,y,z))
    enddo; enddo; enddo
  enddo
#endif
  return
End Subroutine mg3d_$1
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! preprocessor command: recursively generate subroutines for coarser grid 
ifelse($1,$2,,`CREATE_MG(eval($1+1),$2)')') dnl
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! preprocessor command: start generation of grid hierarchy;
!                       finest grid at level 0, coarsest grid at level 9
CREATE_MG(0,12)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! proceed with quasi 2D algorithm when at coarsest level in z
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! preprocessor command: generate subgrid at level xx, z at coarsest level (nz=1)
define(CREATE_MGXY,`
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                          3D-multigrid in xy/ level $1
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Subroutine mg3dxy_$1(f,b)
! call if nz = 1 at current grid level
  Use pe_env
  Use mgctrl
  Implicit None

  Real(kind=drk), Dimension(0:nx/2**$1+1,0:ny/2**$1+1,0:2) :: f,b
  Real(kind=drk), Save, Dimension(0:nx/2**$1+1,0:ny/2**$1+1,0:2) :: g
  Real(kind=drk), Save, &
     Dimension(0:nx/2**($1+1)+1,0:ny/2**($1+1)+1,0:2) :: fc,bc
  Integer :: x,y,z,i,j,nu4,levelz
  Real(kind=drk), Parameter :: dx=lx/(nx*pex)*2**$1, dy=ly/(ny*pey)*2**$1, &
      dz=lz/pez
  Real(kind=drk) :: dtau

! calculate time step dtau for smoother
! cpp commands added to deal with possibility that dz is arbitrary in 2D
! simulation
#if (nz == 1 && pez == 1)
  dtau=1./(c2*(4./dx**2+4./dy**2)+c1)
#else
  dtau=1./(c2*(4./dx**2+4./dy**2+4./dz**2)+c1)
#endif

! calculate z level
  levelz=0
  do while (mod(nz/2**levelz,2) == 0) 
    levelz=levelz+1
  enddo
! pre-smoothing (nu1 iterationsi, choose JACOBI, GS_LEX, GS_RB)
#ifdef JACOBI
  do i=1,nu1
    call bound_phi($1,levelz,f)
    do y=1,ny/2**$1; do x=1,nx/2**$1
      g(x,y,1)=c2*DELTA2(f)-c1*f(x,y,1)+c3*b(x,y,1)
    enddo; enddo
    do y=1,ny/2**$1; do x=1,nx/2**$1
      f(x,y,1)=f(x,y,1)+dtau*g(x,y,1)
    enddo; enddo
  enddo
#endif
#ifdef GS_LEX
  do i=1,nu1
    call bound_phi($1,levelz,f)
    do y=1,ny/2**$1; do x=1,nx/2**$1
      f(x,y,1)=f(x,y,1)+dtau*(c2*DELTA2(f)-c1*f(x,y,1)+c3*b(x,y,1))
    enddo; enddo
  enddo
#endif
#ifdef GS_RB
  do i=1,2*nu1
    call bound_phi($1,levelz,f)
    do y=1,ny/2**$1; do x=1,nx/2**$1
      if (mod(x+my_pex*nx/2**$1+y+my_pey*ny/2**$1,2)==mod(i,2)) &
      	f(x,y,1)=f(x,y,1)+dtau*(c2*DELTA2(f)-c1*f(x,y,1)+c3*b(x,y,1))
    enddo; enddo
  enddo
#endif
! multigrid recursion
  if (mod(nx/2**$1,2)+mod(ny/2**$1,2) == 0) then
! restriction of defect (9 or 4 point stencil), guess for coarse grid solution
    call bound_phi($1,levelz,f)
    do y=1,ny/2**$1; do x=1,nx/2**$1
      g(x,y,1)=-c3*(c2*DELTA2(f)-c1*f(x,y,1)+c3*b(x,y,1))
    enddo; enddo
#ifdef alt_restrict_prolong
    call boundary($1,levelz,g)
    do y=1,ny/2**($1+1); do x=1,nx/2**($1+1)
      bc(x,y,1)=(g(2*x  ,2*y  ,1)   +&
                 g(2*x-1,2*y  ,1)/2 + g(2*x+1,2*y  ,1)/2 +&
                 g(2*x  ,2*y-1,1)/2 + g(2*x  ,2*y+1,1)/2 +&
                 g(2*x-1,2*y-1,1)/4 + g(2*x-1,2*y+1,1)/4 +&
                 g(2*x+1,2*y-1,1)/4 + g(2*x+1,2*y+1,1)/4    )/4
      fc(x,y,1)=0.
    enddo; enddo
#else
    do y=1,ny/2**($1+1); do x=1,nx/2**($1+1)
      bc(x,y,1)=(g(2*x  ,2*y  ,1) + g(2*x-1,2*y  ,1) +&
                 g(2*x  ,2*y-1,1) + g(2*x-1,2*y-1,1))/4
      fc(x,y,1)=0.
    enddo; enddo
#endif
! multigrid recursion
ifelse($1,$2,,`    call mg3dxy_`'eval($1+1)(fc,bc)')
#ifdef WCYCLE
ifelse($1,$2,,`    call mg3dxy_`'eval($1+1)(fc,bc)')
#endif
! subtract prolongation (9 or 4 point stencil)
#ifdef alt_restrict_prolong
    call boundary($1+1,levelz,fc)
    do y=1,ny/2**($1+1); do x=1,nx/2**($1+1)
      f(2*x  ,2*y  ,1) = f(2*x  ,2*y  ,1) - fc(x,y,1)
      f(2*x-1,2*y  ,1) = f(2*x-1,2*y  ,1) - (fc(x,y,1)+fc(x-1,y,1))/2
      f(2*x  ,2*y-1,1) = f(2*x  ,2*y-1,1) - (fc(x,y,1)+fc(x,y-1,1))/2
      f(2*x-1,2*y-1,1) = f(2*x-1,2*y-1,1) &
           -(fc(x,y,1)+fc(x,y-1,1)+fc(x-1,y,1)+fc(x-1,y-1,1))/4
    enddo; enddo
#else
    do y=1,ny/2**($1+1); do x=1,nx/2**($1+1)
      f(2*x  ,2*y  ,1) = f(2*x  ,2*y  ,1) - fc(x,y,1)
      f(2*x-1,2*y  ,1) = f(2*x-1,2*y  ,1) - fc(x,y,1)
      f(2*x  ,2*y-1,1) = f(2*x  ,2*y-1,1) - fc(x,y,1)
      f(2*x-1,2*y-1,1) = f(2*x-1,2*y-1,1) - fc(x,y,1)
    enddo; enddo
#endif
    nu4=nu2
  else
! additional smoothing iterations at coarsest grid
    nu4=nu2+nu3
    if (use_sor) then 
! successive overrelaxation at coarsest grid
      call sor($1,levelz,f,b)
      nu4=0
    endif
  endif
! post-smoothing (nu2 or nu2+nu3 iterations, choose JACOBI, GS_LEX, GS_RB)
#ifdef JACOBI
  do i=1,nu4
    call bound_phi($1,levelz,f)
    do y=1,ny/2**$1; do x=1,nx/2**$1
      g(x,y,1)=c2*DELTA2(f)-c1*f(x,y,1)+c3*b(x,y,1)
    enddo; enddo
    do y=1,ny/2**$1; do x=1,nx/2**$1
      f(x,y,1)=f(x,y,1)+dtau*g(x,y,1)
    enddo; enddo
  enddo
#endif
#ifdef GS_LEX
  do i=1,nu4
    call bound_phi($1,levelz,f)
    do y=1,ny/2**$1; do x=1,nx/2**$1
      f(x,y,1)=f(x,y,1)+dtau*( c2*DELTA2(f)-c1*f(x,y,1)+c3*b(x,y,1) )
    enddo; enddo
  enddo
#endif
#ifdef GS_RB
  do i=1,2*nu4
    call bound_phi($1,levelz,f)
    do y=1,ny/2**$1; do x=1,nx/2**$1
      if (mod(x+my_pex*nx/2**$1+y+my_pey*ny/2**$1,2)==mod(i,2)) &
      	f(x,y,1)=f(x,y,1)+dtau*(c2*DELTA2(f)-c1*f(x,y,1)+c3*b(x,y,1))
    enddo; enddo
  enddo
#endif
  return
End Subroutine mg3dxy_$1
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! preprocessor command: recursively generate subroutines for coarser grid 
ifelse($1,$2,,`CREATE_MGXY(eval($1+1),$2)')') dnl
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! preprocessor command: start generation of grid hierarchy;
!                       finest grid at level 0, coarsest grid at level 9
CREATE_MGXY(0,12)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!------------------------------------------------------------------------------
!                          3D-multigrid / residual
!------------------------------------------------------------------------------

Function res3d(f,b)
! calculate the residual as a measure of the error in the approximate solution

  Use pe_env
  Use mgctrl
  Implicit None
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: f,b
  Real(kind=drk) :: res3d, sum
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), &
       dz=lz/(nz*pez)
  Integer :: x,y,z,pe

  sum=0.
  call bound_phi(0,0,f)
  do z=1,nz; do y=1,ny; do x=1,nx
    sum=sum+abs(c2*DELTA(f)-c1*f(x,y,z)+c3*b(x,y,z))
  enddo; enddo; enddo
   res3d = totalsum(sum)/normalize
  return
End Function res3d

!------------------------------------------------------------------------------
!                        Successive overrelaxation
!------------------------------------------------------------------------------

Subroutine sor(levelxy,levelz,f,b)
  Use pe_env
  Implicit None

  Integer :: x,y,z,n,levelxy,levelz,nnx,nny,nnz
  Real(kind=drk),Dimension(0:nx/2**levelxy+1,0:ny/2**levelxy+1,0:nz/2**levelz+1) :: f,b
  Real(kind=drk), Parameter :: pi2=6.2831853
  Real(kind=drk) :: rjac, omega, dx, dy, dz, cc

! some constants
  dx=lx/(nx*pex)*2**levelxy
  dy=ly/(ny*pey)*2**levelxy
  dz=lz/(nz*pez)*2**levelz

! cpp commands added to deal with possibility that dz is arbitrary in 2D
! simulation
#if (nz == 1 && pez == 1)
  cc=-2*(1/dx**2+1/dy**2)
#else
  cc=-2*(1/dx**2+1/dy**2+1/dz**2)
#endif

! calculate spectral radius of Jacobi iteration (periodic boundary conditions)
#if (nz == 1 && pez == 1) 
  rjac=(cos(pi2*dx/lx)/dx**2+cos(pi2*dy/ly)/dy**2) / (1/dx**2+1/dy**2)
#else
  rjac=(cos(pi2*dx/lx)/dx**2+cos(pi2*dy/ly)/dy**2+cos(pi2*dz/lz)/dz**2) / &
       (1/dx**2+1/dy**2+1/dz**2)
#endif

! calculate overrelaxation parameter
  omega=2/(1+sqrt(1-rjac**2))

! SOR iterations
  do n=1,nu3
    call bound_phi(levelxy,levelz,f)
    do z=1,nz/2**levelz; do y=1,ny/2**levelxy; do x=1,nx/2**levelxy
      nnx=x+my_pex*nx/2**levelxy
      nny=y+my_pey*ny/2**levelxy
      nnz=z+my_pez*nz/2**levelz
      if (mod(nnx+nny+nnz,2)==mod(n,2)) &
                  f(x,y,z)=f(x,y,z)-omega*(DELTA(f)-b(x,y,z))/cc
    enddo; enddo; enddo
  enddo
End Subroutine sor
