!******************************************************************************
!             3D particle code: initialization - support routines
!                            Andreas Zeiler, 1999
!
!                        LATEST CHANGE: July 16, 2001
!
!******************************************************************************

#include "param"

#define nsteps 1000

!------------------------------------------------------------------------------
!                       Maxwellian velocity distribution
!------------------------------------------------------------------------------

Subroutine maxwellian(t,m,vx,vy,vz)
! returns random velocity vector for particle with mass m
! series of random numbers consistent with relativistic / non-relativistic
! Maxwell distribution at temperature t

  Implicit None
  Real(kind=drk) :: t,m,v,f,ekin,gamma,tmp,rx,ry,rz,r,gaussian,rand_num,beta
  Real(kind=prk) :: vx,vy,vz
  Integer :: i
  Logical :: reject
  Logical, Save :: normalize=.true.
  Real(kind=drk), Save :: norm,t_old=0.,m_old=0.

#ifdef relativistic
! relativistic distribution
  gamma(beta)=1./sqrt(1-beta**2)
  ekin(beta)=m*c_2*(gamma(beta)-1)

! Corrected version of the velocity dependence of the relativistic
! Maxwellian.  The whole shebang is f(gamma) = x/K_2(x)*gamma*
! sqrt(gamma^2-1)*exp(-gamma*x) where x is mc^2/kT and K_2 is an
! modified Bessel function.  In the limits of x large, K_2(x) is 
! exp(-x)*sqrt(pi/(2*x)).  Converting to f(v) gives a 
! gamma^5*beta^2.  The result is f(v) = 
! sqrt(2/pi)*(m/t)^1.5*gamma(v)^5*v^2*exp(-ekin(v)/t)

! Since beta varies between 0 and 1, use it instead.  f(beta) = 
! sqrt(2/pi)*(m*c_2/t)^1.5*beta^2/(1-beta^2)^2.5*
! exp(-m*c_2/t*(1/sqrt(1-b^2)-1)).
f(beta) = norm*beta**2*gamma(beta)**5*exp(-ekin(beta)/t)

! First, if m*c^2/t is bigger than 100, calculate a non-relativistic 
! distribution.  Do this to eliminate loss-of precision in large x
! limit.  As a sense of the degree of approximation, the most probable
! velocity for mc^2/t = 100 is 0.14*c.
  if (m*c_2/t .ge. 100.) then
    vx=gaussian()*sqrt(t/m)
    vy=gaussian()*sqrt(t/m)
    vz=gaussian()*sqrt(t/m)
  else 
! normalize distribution such that max(f(beta)) = 1.
    if (m .ne. m_old .or. t .ne. t_old) normalize=.true.
    if (normalize) then
      norm=1.
      tmp=0.
      do i=0,nsteps-1
        tmp=max(tmp,f(real(1.,drk)*i/nsteps))
      enddo
      norm=1./tmp
      t_old=t
      m_old=m
      normalize=.false.
    endif
! randomly generate beta
    reject =.true.
    do while (reject)
      beta=rand_num()
      reject = (rand_num()>f(beta))
    enddo
! randomly generate direction for v
    reject =.true.
    do while (reject)
      rx=rand_num()*2.-1.
      ry=rand_num()*2.-1.
      rz=rand_num()*2.-1.
      r=sqrt(rx**2+ry**2+rz**2)
      reject = (r .ge. 1 .or. r .eq. 0)
    enddo
    v = beta*sqrt(c_2)
    vx=v*rx/r
    vy=v*ry/r
    vz=v*rz/r
  endif
#else
! non-relativistic distribution
  vx=gaussian()*sqrt(t/m)
  vy=gaussian()*sqrt(t/m)
  vz=gaussian()*sqrt(t/m)
#endif
  return
End Subroutine maxwellian

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Subroutine trimaxwellian(t1,t2,t3,m,vx,vy,vz)
! returns random velocity vector for particle with mass m but different temperatures
! along x,y,z. It can be used to get bi-maxwellian distributions with mean field along
! any of the three directions. Right now I've taken out the relativistic part. Will
! add it later if we decide to do relativistic simulations.

  Implicit None
  Real(kind=drk) :: t1,t2,t3,m,gaussian
  Real(kind=prk) :: vx,vy,vz
! non-relativistic distribution
  vx=gaussian()*sqrt(t1/m)
  vy=gaussian()*sqrt(t2/m)
  vz=gaussian()*sqrt(t3/m)
  return

End Subroutine trimaxwellian

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Function gaussian()
! returns a normally distributed deviate with zero mean and unit variance
! (algorithm gasdev from Numerical Recipes)

  Implicit none
  Real(kind=drk) gaussian, rand_num
  Logical, Save :: iset=.true.
  Real(kind=drk), Save :: gset
  Real(kind=drk) :: fac,rsq,v1,v2
  if (iset) then
    rsq=0.
    do while (rsq.ge.1 .or. rsq == 0.)
      v1=2.*rand_num()-1.
      v2=2.*rand_num()-1.
      rsq=v1**2+v2**2
    enddo
    fac=sqrt(-2.*log(rsq)/rsq)
    gset=v1*fac
    gaussian=v2*fac
    iset=.false.
  else
    gaussian=gset
    iset=.true.
  endif
  return
End Function gaussian

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Function rand_num()
! returns a random number
! a single call to the intrinsic "random_number" is incredibly slow on the IBM,
! so we do it more efficiently here

! M. Shay, 10/26/2007: All of the calls to rand_num in p3d expect a 
!                      double precision number, but ran1 returns a single 
!                      precision number. So, make rand_num double precision by
!                      hand. 
!                      WARNING: rand_num not really double precision. 

  Use pe_env
  Use nrtype
  Use nr, Only: ran1
  Implicit none
  Real(kind=drk) rand_num
  Integer, Parameter :: lran=10000
  Real(kind=prk), Save, Dimension(lran) :: rnumbers
  Integer, Save :: num = lran

  if (num>=lran) then
    call ran1(rnumbers)
    num=0
  endif

  num = num+1
  rand_num = real(rnumbers(num),drk)

End Function rand_num
!------------------------------------------------------------------------------
!                       Random particle location
!------------------------------------------------------------------------------

Subroutine location(x,y,z)
! returns random location in space

  use pe_env
  Implicit None
  Real(kind=prk) :: x,y,z
  Real(kind=drk) rand_num

! x=rand_num()*lx; y=rand_num()*ly; z=rand_num()*lz
! x=mod(x,real(lx,prk)); y=mod(y,real(ly,prk)); z=mod(z,real(lz,prk))
  x = lx/pex*(rand_num()+my_pex)
  y = ly/pey*(rand_num()+my_pey)
  z = lz/pez*(rand_num()+my_pez)
End Subroutine location

!------------------------------------------------------------------------------
!                     summation of velocities
!------------------------------------------------------------------------------

Subroutine veloadd(ux,uy,uz,vx,vy,vz)
! add velocity v relativistically / nonrelativistically  correct to vector u

  Use pe_env

  Implicit none
  Real(kind=prk) :: ux,uy,uz
  Real(kind=drk) :: vx,vy,vz,uxp,uyp,uzp,gamma,tmp
  Logical,Save :: message=.true.

#ifdef relativistic
  if (vx**2+vy**2+vz**2 == 0) return
! check whether v is smaller than c, normalize otherwise
  if ((vx**2+vy**2+vz**2) .ge. c_2*0.99) then
    tmp=sqrt(c_2*0.99/(vx**2+vy**2+vz**2))
    vx=vx*tmp; vy=vy*tmp; vz=vz*tmp
    if (message) then
      write(6,*) 'PE ',myproc,' *** warning: velocity exceeds speed of light' 
      call flush(6)
      message=.false.;
    endif
  endif
! calculate parallel component of u
  tmp=(ux*vx+uy*vy+uz*vz)/(vx**2+vy**2+vz**2)
  uxp=tmp*vx; uyp=tmp*vy; uzp=tmp*vz
! calculate inverse of gamma
  gamma=sqrt(1-(vx**2+vy**2+vz**2)/c_2)
! calculate new velocity
  tmp=1+(ux*vx+uy*vy+uz*vz)/c_2
  ux=(ux*gamma+vx+uxp*(1-gamma))/tmp
  uy=(uy*gamma+vy+uyp*(1-gamma))/tmp
  uz=(uz*gamma+vz+uzp*(1-gamma))/tmp
#else
  ux=ux+vx; uy=uy+vy; uz=uz+vz
#endif
  return
End Subroutine veloadd
