!******************************************************************************
!      3D particle code: two stream initialization (homogeneous beams)
!                           Andreas Zeiler, 1999
!
!                       LATEST CHANGE: July 16, 2001
!
!******************************************************************************

#include "param"
#ifndef T_i
#define T_i 0.
#endif
#ifndef T_e
#define T_e 0.
#endif
#ifndef v_0
#define v_0 0.
#endif

Subroutine init_twostream()
  Use pe_env
  Use partfield
  Implicit None

  Real(kind=drk) gaussian, rand_num
  Integer i,x,y,z

  call init_pe_env()
  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** init: compiled for different number of PEs *****'
    call exitallpes()
  endif

! write parameters to log-file
  if (myproc==0) then
     write(6,*) '******** inittwostream *********'
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
     write(6,*) '     part./gridp. = ',ppg
     write(6,*) '     T_i = ',T_i
     write(6,*) '     T_e = ',T_e
     write(6,*) '     v_0 = ',v_0
#ifdef relativistic
     write(6,*) '     relativistic initialization'
#endif
     write(6,*) '********** parameters **********'
  endif

! set electric and magnetic field
  do z=1,nz; do y=1,ny; do x=1,nx
    b1x(x,y,z)=0.; b1y(x,y,z)=0.; b1z(x,y,z)=0.
    e1x(x,y,z)=0.; e1y(x,y,z)=0.; e1z(x,y,z)=0.
  enddo; enddo; enddo 
  np_i=0; np_e=0
! load particles
  do i=1,ppg*nx*ny*nz
    np_i=np_i+1
    if (np_i .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
    endif
    call location(rvi(1,np_i),rvi(2,np_i),rvi(3,np_i))
    call maxwellian(real(T_i,drk), real(1.,drk), &
                    rvi(4,np_i), rvi(5,np_i), rvi(6,np_i))
    call veloadd(rvi(4,np_i),rvi(5,np_i),rvi(6,np_i), &
                 real( v_0,drk), real(0.,drk), real(0.,drk))
  enddo
  do i=1,ppg*nx*ny*nz
    np_e=np_e+1
    if (np_e .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
    endif
    call location(rve(1,np_e),rve(2,np_e),rve(3,np_e))
    call maxwellian(real(T_e,drk), real(m_e,drk), &
                    rve(4,np_e), rve(5,np_e), rve(6,np_e))
    call veloadd(rve(4,np_e),rve(5,np_e),rve(6,np_e), &
                 real(-v_0,drk), real(0.,drk), real(0.,drk))
  enddo
  n_avg=ppg
  call redistribute(rvi,np_i)
  call redistribute(rve,np_e)
  call output(0,.true.)
End Subroutine init_twostream
