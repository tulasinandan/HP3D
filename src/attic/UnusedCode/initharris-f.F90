!******************************************************************************
!    3D particle code: double harris sheet initialization, two-fluid version
!                           Andreas Zeiler, 2001
!
!                        LATEST CHANGE: July 17, 2001
!
!******************************************************************************

#include "param"
#ifndef b0
#define b0 0.
#endif
#ifndef psi0
#define psi0 0.
#endif
#ifndef w0
#define w0 1.
#endif
#ifndef T_i
#define T_i 0.
#endif
#ifndef T_e
#define T_e 0.
#endif
#ifndef n_0
#define n_0 0.
#endif

Subroutine init_harris()
  Use pe_env
  Use partfield
  Implicit None

  Real(kind=drk) gx,gy,gz,yy,b,n,gaussian,sech
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez),&
                               ppi2=6.2831853, ppi4=2*ppi2
  Integer i,x,y,z,pe,sign

  gx(x)=((x-0.5)+my_pex*nx)*dx
  gy(y)=((y-0.5)+my_pey*ny)*dy
  gz(z)=((z-0.5)+my_pez*nz)*dz
  sech(yy)=1/cosh(yy)

! shape of the sheet
  b(yy)=tanh((yy-ly*.25)/w0)-tanh((yy-ly*.75)/w0)-1
  n(yy)=(sech(real(((yy-ly*.25)/w0),drk))**2+ &
         sech(real(((yy-ly*.75)/w0),drk))**2)/(2*(T_i+T_e))

  call init_pe_env()
  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** init-hybrid: compiled for different number of PEs *****'
    call exitallpes()
  endif

! write parameters to log-file
  if (myproc==0) then
    write(6,*) '********** initharris **********'
    write(6,*) '***** (two-fluid version) ******'
    write(6,*) '********** parameters **********'
    write(6,*) '     grid points nx = ',nx
    write(6,*) '     grid points ny = ',ny
    write(6,*) '     grid points nz = ',nz
    write(6,*) '     processors pex = ',pex
    write(6,*) '     processors pey = ',pey
    write(6,*) '     processors pez = ',pez
    write(6,*) '     edge length lx = ',lx
    write(6,*) '     edge length ly = ',ly
    write(6,*) '     edge length lz = ',lz
    write(6,*) '     T_i = ',T_i
    write(6,*) '     T_e = ',T_e
    write(6,*) '     n_0 = ',n_0
    write(6,*) '     b0 = ',b0
    write(6,*) '     w0 = ',w0
    write(6,*) '     psi0 = ',psi0
    write(6,*) '********** parameters **********'
  endif

! set magnetic field, density, ion and electron pressure and ion current
  do z=1,nz; do y=1,ny; do x=1,nx
    b1x(x,y,z)=b(gy(y)); b1y(x,y,z)=0.; b1z(x,y,z)=b0
    n1(x,y,z)=n_0+n(gy(y)) 
    pe1(x,y,z)=n1(x,y,z)*T_e; pi1(x,y,z)=n1(x,y,z)*T_i
    sign=1; if (gy(y) .gt. ly/2) sign=-1
    j1x(x,y,z)=0.; j1y(x,y,z)=0.; j1z(x,y,z)=-(n(gy(y)))*2*T_i/w0*sign
  enddo; enddo; enddo
! add perturbation to magnetic field to start reconnection
  do z=1,nz; do y=1,ny; do x=1,nx
    b1x(x,y,z)=b1x(x,y,z)-psi0*ppi4/ly*sin(ppi2*gx(x)/lx)*sin(ppi4*gy(y)/ly)
    b1y(x,y,z)=b1y(x,y,z)+psi0*ppi2/lx*cos(ppi2*gx(x)/lx)*(1-cos(ppi4*gy(y)/ly))
  enddo; enddo; enddo
! write start-up file
  call output(0,.true.)
End Subroutine init_harris
