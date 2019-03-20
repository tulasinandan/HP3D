!******************************************************************************
!             3D particle code: Ion-Acoustic Wave setup
!                           template by Andreas Zeiler, 1999
!                           Modified for wave by Michael Shay, 2004
!
!                      LATEST CHANGE: June 24, 2004
! M. Shay: Early July, 2004: Wrote ion-acoustic wave setup.

! M. Shay: 7/29/2004:        Added debeye length corrections to ion
!                            acoustic wave setup. 
!
! M. Shay: 11/5/2004:        Added a constant velocity to be added to
!                            electrons and ions so that simulation is
!                            done in the wave frame. 
!******************************************************************************



! initperturb: current along z
! initperturb2: current along x

#include "param"
#define initperturb 1
#define initperturb2 2
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
#define n_0 1.
#define no_background_density
#endif
#ifndef v_waveadd 
#define v_waveadd 0.
#endif

Subroutine init_waveE()
  Use pe_env
  Use partfield
  Implicit None

  Real(kind=drk) gx,gy,gz,yy,b,rand_num,sech,nsup,j,ntwid,c_s,vx,xx
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez),&
                          pi2=6.2831853, pi4=2.*pi2
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: curl_x,curl_y,curl_z,gridn
  Integer i,x,y,z,pe
  Real(kind=drk) :: realx,realy,realz
  Integer :: xval,yval,zval,ierror
  Integer :: numpart,numparti,numparte,numpartgrid

  Real(kind=drk) :: xp,yp,zp,xm,ym,zm,inter_n,inter_jx,inter_jy,inter_jz
  Real(kind=drk) :: maxn,sumn,globalmaxn,globalsumn
  Real(kind=drk) :: edebye,debcorr,kval,ni,ne,vix,vex,Eval
  Real(kind=drk) :: edebye2,idebye2,kval2,c_se2,c_si2,omega,omega2
  Real(kind=drk) :: const1,const2,nemag
  Real(kind=drk) :: sumni,sumne,globalsumni,globalsumne

! print *
!!!real(kind=drk) :: valtemp(nx)



! Runs from 0 to lx basically over the whole simulation.
  gx(x)=((x-0.5)+my_pex*nx)*dx
  gy(y)=((y-0.5)+my_pey*ny)*dy
  gz(z)=((z-0.5)+my_pez*nz)*dz
  sech(yy)=1/cosh(yy)

! 1D equilibrium
! Don't worry. c_si2 and other variables are defined after vix(xx), but vix(xx) is a little
! function which is not executed until after c_si2 is defined. 

  ntwid(xx) = npert*sin(kval*xx)
  ni(xx) = n_0 + ntwid(xx)
  ne(xx) = n_0 + ntwid(xx)*nemag
  vix(xx) = ntwid(xx)*omega/kval + (v_waveadd)
  vex(xx) = ntwid(xx)*nemag*omega/kval + (v_waveadd)
  Eval(xx) = (1.*npert/n_0)*cos(kval*xx)* &
             (-omega2/kval + kval*c_si2)

! Some variables for the equilibrium

  edebye2 = 1.*T_e/c_2 ; idebye2 = 1.*T_i/c_2
  kval = 6.28318530718/lx ; kval2 = kval**2
  c_se2 = 1.*T_e*gamma_sounde ; c_si2 = 1.*T_i*gamma_soundi
  const1 = kval2*c_se2*(1 + c_si2/c_se2 + kval2*idebye2)
  const2 = 1. + kval2*edebye2 + m_e + m_e*kval2*idebye2
  omega = sqrt(const1/const2) ; omega2 = const1/const2
  nemag = 1. - omega2/c_2 + kval2*idebye2

!  b(yy)=tanh((yy-ly*.25)/w0)-tanh((yy-ly*.75)/w0)+&
!	tanh((yy-ly*1.25)/w0)-tanh((yy+ly*.25)/w0)+1
!  j(yy) = -(sech((yy-ly*.25)/w0)**2-sech((yy-ly*.75)/w0)**2+&
!        sech((yy-ly*1.25)/w0)**2-sech((yy+ly*.25)/w0)**2)/w0
!  n(yy) = (1-(tanh((yy-ly*.25)/w0)-tanh((yy-ly*.75)/w0)+&
!      tanh((yy-ly*1.25)/w0)-tanh((yy+ly*.25)/w0)+1)**2)/(2*(T_i+T_e))
  nsup = 1.

! Anything processor dependent must be done after this next line.
  call init_pe_env()
  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** init: compiled for different number of PEs *****'
    call exitallpes()
  endif

!!!  do x=1,nx
!!!    valtemp(x) = vix(gx(x))
!!!  end do

!!!  print *,'myproc = ',myproc,'  vix = ',valtemp

  ! For the purposes of determining how much density each particle
  ! represents, determine maximum density.
  ! For initializing particles, determine maximum density.

       ! Note: Becuase this is a linear wave, even though ni(xx) is
       ! not the same as ne(xx), their sum over all space is the
       ! same. 
  sumn = 0. ; globalsumn = 0. ; maxn = 0. ; globalmaxn = 0.
  do z=1,nz ; do y=1,ny ; do x=1,nx
    sumni = sumni + ni(gx(x))
    sumne = sumne + ne(gx(x))
!    maxn = max(maxn,ni(gx(x)))
  enddo ; enddo ; enddo

  globalsumni=totalsum(sumni)
  globalsumne=totalsum(sumne)

  numpart=ppg*nx*pex*ny*pey*nz*pez

! write parameters to log-file
  if (myproc==0) then
#if (init_scheme == initperturb)
    write(6,*) '********** perturb **********'
#elif (init_scheme == initperturb2)
    write(6,*) '********** initperturb2 *********'
#endif
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
    write(6,*) '     n_0 = ',n_0
    write(6,*) '     npert = ',npert
    write(6,*) '     gamma_soundi = ',gamma_soundi
    write(6,*) '     gamma_sounde = ',gamma_sounde
#ifdef relativistic
    write(6,*) '     relativistic initialization'
#endif
    write(6,*) '********** parameters **********'
  endif

! set electric and magnetic field (no kinetic equilibrium with guide field bz)
  do z=1,nz; do y=1,ny; do x=1,nx
    b1x(x,y,z)=0.; b1y(x,y,z)=0. ; b1z(x,y,z)=0.
    e1x(x,y,z)=Eval(gx(x)); e1y(x,y,z)=0.; e1z(x,y,z)=0.
  enddo; enddo; enddo 


!this updates the ghost_cells x=0,x=nx+1,...
  call bound_b(b1x,b1y,b1z)
  do z=1,nz; do y=1,ny; do x=1,nx
    curl_x(x,y,z)=(b1z(x,y+1,z)-b1z(x,y-1,z))/(2.*dy)- &
                  (b1y(x,y,z+1)-b1y(x,y,z-1))/(2.*dz) 
    curl_y(x,y,z)=(b1x(x,y,z+1)-b1x(x,y,z-1))/(2.*dz)- &
                  (b1z(x+1,y,z)-b1z(x-1,y,z))/(2.*dx) 
    curl_z(x,y,z)=(b1y(x+1,y,z)-b1y(x-1,y,z))/(2.*dx)- &
                  (b1x(x,y+1,z)-b1x(x,y-1,z))/(2.*dy)
  enddo; enddo; enddo
  call bound_j(curl_x,curl_y,curl_z)

  np_i=0; np_e=0

! Load particles. This is an attempt to do a relatively quiet start.
  ! The approximate number of particles necessary to give the wanted
  ! density are loaded uniformly around a particular grid point. This
  ! process is repeated for all grid points. 

! load ions
do z=1,nz ; do y=1,ny ; do x=1,nx
  numpartgrid=nint(ni(gx(x))/globalsumni*numpart)
  do i=1,numpartgrid
     np_i=np_i+1
     if (np_i .gt. maxparticles) then
       write(6,*) '***** init: particle buffer overflow *****'
       call exitallpes()
     endif
     rvi(1,np_i)=gx(x) + (rand_num()-0.5)*dx
     rvi(2,np_i)=gy(y) + (rand_num()-0.5)*dy
     rvi(3,np_i)=gz(z) + (rand_num()-0.5)*dz
     call maxwellian(real(T_i,drk), real(1.,drk), &
                      rvi(4,np_i), rvi(5,np_i), rvi(6,np_i))
     call veloadd(rvi(4,np_i),rvi(5,np_i),rvi(6,np_i), &
        vix(real(rvi(1,np_i),drk)),real(0.,drk),real(0.,drk) )
  end do
end do ; end do ; end do

! load electrons
do z=1,nz ; do y=1,ny ; do x=1,nx
  numpartgrid=nint(ne(gx(x))/globalsumne*numpart)
  do i=1,numpartgrid
     np_e=np_e+1
     if (np_e .gt. maxparticles) then
       write(6,*) '***** init: particle buffer overflow *****'
       call exitallpes()
     endif
     rve(1,np_e)=gx(x) + (rand_num()-0.5)*dx
     rve(2,np_e)=gy(y) + (rand_num()-0.5)*dy
     rve(3,np_e)=gz(z) + (rand_num()-0.5)*dz
     call maxwellian(real(T_e,drk), real(m_e,drk), &
                      rve(4,np_e), rve(5,np_e), rve(6,np_e))
     call veloadd(rve(4,np_e),rve(5,np_e),rve(6,np_e), &
        vex(real(rve(1,np_e),drk)),real(0.,drk),real(0.,drk) )
  end do
end do ; end do ; end do

! Calculate the number of particles per density. There is currently
! a small accuracy error that occurs if the total number of electrons
! does not equal the total number of ions. 

  numparti=totalsumint(np_i)
  numparte=totalsumint(np_e)

  print *,myproc,'   np_i = ',np_i,'   np_e = ',np_e

  print *,myproc,'  numparti = ',numparti,'   numparte = ',numparte

  n_avg = 1.*numparti/globalsumni

  print *,myproc,'   n_avgi = ',1.*numparti/globalsumni,'   n_avge = ',  &
          1.*numparte/globalsumne

!#ifndef subtract_average_rho
!! force global charge neutrality
!  np_i=min(np_i,np_e); np_e=min(np_i,np_e)
!print *,'np_i = ',np_i
!print *,'np_e = ',np_e
!#endif

  call redistribute(rvi,np_i)
  call redistribute(rve,np_e)

  call output(0,.true.)
End Subroutine init_waveE

