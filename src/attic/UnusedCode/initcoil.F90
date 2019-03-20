!******************************************************************************
!        3D particle code: Dave initialization
!                           Andreas Zeiler, 2000
!
!                       LATEST CHANGE: July 17, 2001
!
!******************************************************************************

! based on initelectron3: current along z and homogeneous ions

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
#ifndef T_e
#define T_e 0.
#endif
#ifndef T_i
#define T_i 0.
#endif

Subroutine init_coil()
  Use pe_env       ! in modules.F90, has PE info
  Use partfield    ! in modules.F90
  Implicit None

  Real(kind=drk) gx,gy,gz,yy,rx,ry,b,n,ntotal,ey_bx,rand_num
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez),&
                               pi2=6.2831853, pi4=2*pi2,ddy=dy*0.01
  Integer i,x,y,z,pe,sign,nn

  gx(x)=((x-0.5)+my_pex*nx)*dx
  gy(y)=((y-0.5)+my_pey*ny)*dy
  gz(z)=((z-0.5)+my_pez*nz)*dz

  n(rx,ry)= exp(-((rx-0.3)/lx)**2-((ry-0.75)/ly)**2)  ! initial density

  call init_pe_env()
  write(6,*) 'into init_dave()'
  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** init: compiled for different number of PEs *****'
    call exitallpes()
  endif

! write parameters to log-file
  if (myproc==0) then
    write(6,*) '******** initdave ********'
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
    write(6,*) '     m_e = ',m_e
    write(6,*) '     c_2 = ',c_2
    write(6,*) '     part./gridp. = ',ppg
    write(6,*) '     T_e = ',T_e
    write(6,*) '     T_i = ',T_i
    write(6,*) '     b0 = ',b0
    write(6,*) '     w0 = ',w0
    write(6,*) '     psi0 = ',psi0
#ifdef relativistic
    write(6,*) '     relativistic initialization'
#endif
    write(6,*) '********** parameters **********'
  endif

! set electric and magnetic field
  do z=1,nz; do y=1,ny; do x=1,nx
    b1x(x,y,z)=0.; b1y(x,y,z)=0.; b1z(x,y,z)=b0
    e1x(x,y,z)=0.; e1y(x,y,z)=0.; e1z(x,y,z)=0.
  enddo; enddo; enddo 

  np_i=0; np_e=0

! load electrons
  do i=1,nx*ny*nz*ppg*20
    np_e=np_e+1
    if (np_e .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
    endif
    call location(rve(1,np_e),rve(2,np_e),rve(3,np_e))
!        accept or reject particle according to local density
!    if (rand_num()*20 .ge. n( real(rve(1,np_e),drk), real(rve(2,np_e),drk)) ) then
    if ( (real(rve(1,np_e),drk)-coilx)**2 + (real(rve(2,np_e),drk)-coily)**2 .ge. coilrad**2 ) then
      np_e=np_e-1
    else
       rve(4,np_e)=0
       rve(5,np_e)=0
       rve(6,np_e)=1    !all particles flowing in z direction
!      call maxwellian(real(T_e,drk), real(m_e,drk), &
!                      rve(4,np_e), rve(5,np_e), rve(6,np_e))
    endif
  enddo

  n_avg=ppg
  call redistribute(rve,np_e)

! load ions
  ntotal=np_e        ! quasi-neutrality
  ntotal=totalsum(ntotal)
  np_i=int(ntotal/(pex*pey*pez)+0.5)
  if (myproc == 0) np_i=int(ntotal-np_i*(pex*pey*pez-1)+0.5)
  if (np_i .gt. maxparticles) then
    write(6,*) '***** init: particle buffer overflow *****'
    call exitallpes()
  endif
  do i=1,np_i
    rvi(1,i) = rve(1,i)
    rvi(2,i) = rve(2,i)
    rvi(3,i) = rve(3,i)
    rvi(4,i)=0
    rvi(5,i)=0
    rvi(6,i)=0
!    call maxwellian(real(T_i,drk), real(1.,drk), rvi(4,i), rvi(5,i), rvi(6,i))
  enddo
  call redistribute(rvi,np_i)

! calculate electric field
  call calc_rho(rvi,np_i,real(+1.,drk),.true.)
  call calc_rho(rve,np_e,real(-1.,drk),.false.)
  call subtract_average(rho)
  call smooth_rho(rho,rho)
  call divergence(e1x,e1y,e1z,b1x,b1y,b1z,.false.)

! write output file
  call output(0,.true.)

End Subroutine init_coil



