!******************************************************************************
!                    3D particle code: consistency check
!                           Andreas Zeiler, 1999
!
!                      LATEST CHANGE: July 13, 2001
!
!******************************************************************************

#include "param"

#define periodic 0
#define gem_reconnection 1
#define gem_3d 2

#define initrandom 0
#define initharris 1
#define initgem 2
#define inittest 3
#define initharris2 4
#define initelectron 5
#define initelectron2 6
#define initelectron3 7
#define inittwostream 8
#define initelectron4 9

Program consistency_check
! check consistency of parameter file
  Real :: tmp1,tmp2,tmp3  
  Integer :: tmp4

#if (init_scheme == initgem)
#if (boundary_condition != gem_reconnection && boundary_condition != gem_3d)
  write(6,*) '*** init routine "initgem" should be combined'
  write(6,*) '    with "boundary_condition gem_reconnection"'
  write(6,*) '    or "boundary_condition gem_3d"'
#endif
#endif

#if (boundary_condition == periodic)
#ifndef alt_restrict_prolong
  write(6,*) '*** for better multigrid performance you might'
  write(6,*) '    want to use the switch "alt_restrict_prolong"'
#endif
#else
#ifdef alt_restrict_prolong
  write(6,*) '*** use switch "alt_restrict_prolong" only in'
  write(6,*) '    combination with "boundary_condition periodic"'
#endif
#ifdef hybrid
  write(6,*) '*** hybrid code currently implemented with'
  write(6,*) '    periodic boundary conditions only'
#endif
#ifdef twofluid
  write(6,*) '*** twofluid code currently implemented with'
  write(6,*) '    periodic boundary conditions only'
#endif
#endif
#if (boundary_condition == gem_reconnection)
#if (nz != 1 || pez != 1)
  write(6,*) '*** boundary conditions "gem_reconnection" are'
  write(6,*) '    inconsistent in three-dimensional case'
  write(6,*) '    for 2-D run set "nz 1" and "pez 1"'
#endif
#endif

#if (init_scheme == initelectron)
#ifndef subtract_average_rho
  write(6,*) '*** init routine "initelectron" must be combined'
  write(6,*) '    with the switch "subtract_average_rho"'
#endif
#endif

#if (init_scheme == initelectron2)
#ifndef subtract_average_rho
  write(6,*) '*** init routine "initelectron2" must be combined'
  write(6,*) '    with the switch "subtract_average_rho"'
#endif
#endif

#if (init_scheme == initelectron3)
#ifndef subtract_average_rho
  write(6,*) '*** init routine "initelectron3" must be combined'
  write(6,*) '    with the switch "subtract_average_rho"'
#endif
#endif

#if (init_scheme == initelectron4)
#ifndef subtract_average_rho
  write(6,*) '*** init routine "initelectron4" must be combined'
  write(6,*) '    with the switch "subtract_average_rho"'
#endif
#endif

#if (init_scheme == inittwostream)
#ifndef electrostatic
  write(6,*) '*** init routine "inittwostream" must be combined'
  write(6,*) '    with the switch "electrostatic"'
#endif
#endif

  tmp1=nx*pex/lx
  tmp2=ny*pey/ly
  tmp3=nz*pez/lz
  if (tmp1/tmp2 .gt. 1.2 .or. tmp2/tmp1 .gt. 1.2) then
    write(6,*) '*** resolutions in x and y differ substantially'
    write(6,*) '    presumably bad multigrid performance'
  endif
  if (nz.ne.1 .or. pez.ne.1) then
    if (tmp1/tmp3 .gt. 1.2 .or. tmp3/tmp1 .gt. 1.2) then
      write(6,*) '*** resolutions in x and z differ substantially'
      write(6,*) '    presumably bad multigrid performance'
    endif
    if (tmp2/tmp3 .gt. 1.2 .or. tmp3/tmp2 .gt. 1.2) then
      write(6,*) '*** resolutions in y and z differ substantially'
      write(6,*) '    presumably bad multigrid performance'
    endif
  endif
 
#ifdef n_0
  if (n_0 .le. 1e-5) then
      write(6,*) '*** to run without background density "n_0"'
      write(6,*) '    must be undefined'
  endif
#endif
  tmp4=nx
  do while(mod(tmp4,2)==0)
    tmp4=tmp4/2
  enddo
  if (tmp4 .ne. 1) then
    write(6,*) '*** for best multigrid performance nx should be a power of 2'
  endif
  tmp4=ny
  do while(mod(tmp4,2)==0)
    tmp4=tmp4/2
  enddo
  if (tmp4 .ne. 1) then
    write(6,*) '*** for best multigrid performance ny should be a power of 2'
  endif
  tmp4=nz
  do while(mod(tmp4,2)==0)
    tmp4=tmp4/2
  enddo
  if (tmp4 .ne. 1) then
    write(6,*) '*** nz must be a power of 2'
  endif
#ifdef particle
#ifndef coarsegrid_sor
  write(6,*) '*** for better multigrid performance you might'
  write(6,*) '    want to use the switch "coarsegrid_sor"'
#endif
#ifndef skip_poisson
  write(6,*) '*** for better performance you might'
  write(6,*) '    want to use the switch "skip_poisson"'
#endif
#endif
#ifndef miniter
  write(6,*) '*** for safer multigrid convergence you might'
  write(6,*) '    want to use the switch "miniter"'
#endif
#ifndef JACOBI
#ifndef GS_LEX
#ifndef GS_RB
  write(6,*) '*** using JACOBI smoother for multigrid (default)'
#endif
#endif
#endif

#ifdef hybrid
#ifdef twofluid
  write(6,*) '*** "twofluid" and "hybrid" must not be used simultaneously'
#endif
#endif

End Program consistency_check
