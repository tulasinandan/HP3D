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
!
! M. Shay: 12/13/2004:       Electric field and thus electron density are
!                            now calculated by taking E = -Te grad(ni)/ni.
!                            This includes nonlinearities in the electric field. 
!
! M. Shay: 3/19/2006:        Error in previous initializations:
!                            gamma_soundi should equal 3. Debye
!                            length should have a gamma in it. Ion
!                            temperature should have a spatial
!                            dependence. Fixed these errors. 
!
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

Subroutine init_waveH()
  Use pe_env
  Use partfield
  Implicit None

  Real(kind=drk) gx,gy,gz,yy,b,rand_num,sech,j,ntwid,c_s,vx,xx
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez),&
                          pi2=6.2831853, pi4=2.*pi2
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: curl_x,curl_y,curl_z,gridn
  Integer i,x,y,z,pe,ivx,ivy,ivz
  Real(kind=drk) :: realx,realy,realz
  Integer :: xval,yval,zval,ierror
  Integer :: numpart,numparti,numparte,numpartgrid,randi

  Real(kind=drk) :: xp,yp,zp,xm,ym,zm,inter_n,inter_jx,inter_jy,inter_jz
  Real(kind=drk) :: maxn,sumn,globalmaxn,globalsumn
  Real(kind=drk) :: edebye,debcorr,kval,ni,ne,vix,vex,Eval,dxni,dxEval
  Real(kind=drk) :: edebye2,idebye2,kval2,c_se2,c_si2,omega,omega2
  Real(kind=drk) :: const1,const2,nemag
  Real(kind=drk) :: sumni,sumne,globalsumni,globalsumne, Pi, Ti
  Real(kind=drk) :: feps=0.00001,sumfvali,sumfvale,globalsumfvali &
                   ,globalsumfvale,dvxi,dvyi,dvzi,dvxe,dvye,dvze &
                   ,residuei,residuee,totalresiduei,totalresiduee &
                   ,minresiduei,minresiduee,maxresiduei,maxresiduee &
                   ,numpartgridflt,fval,fvaly,fvalz
  Real(kind=drk) :: vymaxi,vzmaxi,vymaxe,vzmaxe
  Real(kind=drk) :: vxvali(nvxi),vyvali(nvyi),vzvali(nvzi) &
                   ,vxvale(nvxe),vyvale(nvye),vzvale(nvze) 





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
  vix(xx) = ntwid(xx)*omega/kval + (v_waveadd)
  dxni(xx) = kval*npert*cos(kval*xx)
  Eval(xx) = -1.*(T_e)*dxni(xx)/ni(xx)
  dxEval(xx) = (T_e)*kval2*npert* &
      (sin(kval*xx)*ni(xx) + npert*(cos(kval*xx))**2)/(ni(xx))**2
  ne(xx) = ni(xx) - dxEval(xx)/c_2
  vex(xx) = vix(xx)
  Pi(xx) = T_i*n_0*(1. + gamma_soundi*ntwid(xx)/n_0) 
  Ti(xx) = Pi(xx)/ni(xx)

! Some variables for the equilibrium

  edebye2 = 1.*gamma_sounde*T_e/c_2 ; idebye2 = 1.*gamma_soundi*T_i/c_2
  kval = 6.28318530718/lx ; kval2 = kval**2
  c_se2 = 1.*T_e*gamma_sounde ; c_si2 = 1.*T_i*gamma_soundi
  const1 = kval2*c_se2*(1 + c_si2/c_se2 + kval2*idebye2)
  const2 = 1. + kval2*edebye2 + m_e + m_e*kval2*idebye2
  omega = sqrt(const1/const2) ; omega2 = const1/const2
  nemag = 1. - omega2/c_2 + kval2*idebye2

! Anything processor dependent must be done after this next line.
  call init_pe_env()
  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** init: compiled for different number of PEs *****'
    call exitallpes()
  endif

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


! Calculate integrated density. Necessary to determine n_avg,
! which is the normalization of the distribution function. 
  sumn = 0. ; globalsumn = 0. ; maxn = 0. ; globalmaxn = 0.
  do z=1,nz ; do y=1,ny ; do x=1,nx
    sumni = sumni + ni(gx(x))
    sumne = sumne + ne(gx(x))
!    maxn = max(maxn,ni(gx(x)))
  enddo ; enddo ; enddo

  globalsumni=totalsum(sumni)
  globalsumne=totalsum(sumne)

! Total number of particles
numpart=ppg*nx*pex*ny*pey*nz*pez
np_i=0; np_e=0

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


! The variables that you need are: 

!    disti  - distribution function in x,vx space
!    Ti     - temperature.

! feps = value of dist function at vymax and vzmax.

! Max values of v
vymaxi = sqrt(-2.*T_i*log(feps))
vzmaxi = vymaxi
vymaxe = sqrt(-2.*T_e/m_e*log(feps))
vzmaxe = vymaxe

print *,'myproc,vymaxi,vzmaxi = ',myproc,vymaxi,vzmaxi

! Grid spacing along x,y,z, vx, vy, vz

! dx,dy,dz determined previously
dvxi = 2.*vxmaxi/(nvxi-1) ; dvxe = 2.*vxmaxe/(nvxe-1)
dvyi = 2.*vymaxi/(nvyi-1) ; dvzi = 2.*vzmaxi/(nvzi-1)
dvye = 2.*vymaxe/(nvye-1) ; dvze = 2.*vzmaxe/(nvze-1)

! Arrays which contain V values
vxvali = (/ (-vxmaxi+(i-1)*dvxi, i=1,nvxi ) /)
vxvale = (/ (-vxmaxe+(i-1)*dvxe, i=1,nvxe ) /)
vyvali = (/ (-vymaxi+(i-1)*dvyi, i=1,nvyi ) /)
vzvali = (/ (-vzmaxi+(i-1)*dvzi, i=1,nvzi ) /)
vyvale = (/ (-vymaxe+(i-1)*dvye, i=1,nvye ) /)
vzvale = (/ (-vzmaxe+(i-1)*dvze, i=1,nvze ) /)

! Zero integral of fval
sumfvali = 0. ; sumfvale = 0.

!############################## Initialize ions
! The Ti^3/2 is necessary to renormalize the velocity integrations
! so that they are independent of x. 
do ivx=1,nvxi ; do x=1,nx
  disti(x,ivx) = ni(gx(x))/(Ti(gx(x)))**(real(1.5,kind=drk)) &
                *exp(-(vxvali(ivx)-vix(gx(x)))**2/(2*Ti(gx(x))/1.)) ! m_i=1
end do ; end do

! Integrate ion dist. to determine maximum value.
! WARNING: The order of summation is very important. The values summed first have lower accuracy. 
! That is why x is summed last.
do x=1,nx ; do y=1,ny ; do z=1,nz ; do ivx=1,nvxi ; do ivy=1,nvyi ; do ivz=1,nvzi
  ! Ion distribution functions along vy and vz.
  fvaly = exp(-(vyvali(ivy))**2/(2*Ti(gx(x))/1.)) ! m_i=1
  fvalz = exp(-(vzvali(ivz))**2/(2*Ti(gx(x))/1.)) !m_i=1

  ! total distribution function
  fval = disti(x,ivx)*fvaly*fvalz

  sumfvali = sumfvali + fval
end do ; end do ; end do ; end do ; end do ; end do

globalsumfvali = totalsum(sumfvali)

residuei = 0.  ! redisdue of extra non-integer particles

! Load ions
do x=1,nx ; do y=1,ny ; do z=1,nz ; do ivx=1,nvxi ; do ivy=1,nvyi ; do ivz=1,nvzi
  ! Ion distribution functions along vy and vz.
  fvaly = exp(-(vyvali(ivy))**2/(2*Ti(gx(x))/1.)) ! m_i=1
  fvalz = exp(-(vzvali(ivz))**2/(2*Ti(gx(x))/1.)) !m_i=1

  ! total distribution function
  fval = disti(x,ivx)*fvaly*fvalz

  numpartgridflt = fval*numpart/globalsumfvali+residuei
  numpartgrid= nint(numpartgridflt)
  residuei = numpartgridflt - numpartgrid

! Put particles on grid
  do i=1,numpartgrid
    np_i=np_i+1
    if (np_i .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
     endif
     rvi(1,np_i)=gx(x) + (rand_num()-0.5)*dx
     rvi(2,np_i)=gy(y) + (rand_num()-0.5)*dy
     rvi(3,np_i)=gz(z) + (rand_num()-0.5)*dz
     rvi(4,np_i)=vxvali(ivx) + (rand_num()-0.5)*dvxi
     rvi(5,np_i)=vyvali(ivy) + (rand_num()-0.5)*dvyi
     rvi(6,np_i)=vzvali(ivz) + (rand_num()-0.5)*dvzi
  end do
end do ; end do ; end do ; end do ; end do ; end do

! Calculate total ion residue

print *,'myproc = ',myproc,'  residuei = ',residuei

totalresiduei = totalsum(residuei)
minresiduei = totalmin(residuei)
maxresiduei = totalmax(residuei)

! Quick and dirty way to get rid of extra residue.
! If the total residue is positive, take the processor with the
! minimum residue and randomly distribute the residue of particles
! over its grid. 
! If the total residue is negative, take the processor with the
! maximum residue and remove the residue of particles from its grid. 

if (totalresiduei >= 0.5 ) then
  ! Add particles to processor with least residue
  if (residuei <= minresiduei) then
    print *,'myproc = ',myproc,'   Adding residue particles'
    do i=1,nint(abs(totalresiduei))
      np_i=np_i+1
      if (np_i .gt. maxparticles) then
        write(6,*) '***** init: particle buffer overflow *****'
        call exitallpes()
      endif
      rvi(1,np_i)=gx(1) + (rand_num())*(gx(nx)-gx(1))
      rvi(2,np_i)=gy(1) + (rand_num())*(gy(ny)-gy(1))
      rvi(3,np_i)=gz(1) + (rand_num())*(gz(nz)-gz(1)) 
      rvi(4,np_i)=-vxmaxi + 2.*rand_num()*vxmaxi
      rvi(5,np_i)=-vymaxi + 2.*rand_num()*vymaxi
      rvi(6,np_i)=-vzmaxi + 2.*rand_num()*vzmaxi
    end do
  end if
else if (totalresiduei <= -0.5) then
  if (residuei >=maxresiduei) then
    print *,'myproc = ',myproc,'   Deleting residue particles'
    do i=1,nint(abs(totalresiduei))
      randi = nint(rand_num()*np_i)
      rvi(:,randi:np_i-1)=rvi(:,randi+1:np_i) ! Shift array
      np_i = np_i - 1
    end do
  end if
else
  if (myproc == 0) print *,'**** No Residue ****'
end if 

! Load particles. This is an attempt to do a relatively quiet start.
  ! The approximate number of particles necessary to give the wanted
  ! density are loaded uniformly around a particular grid point. This
  ! process is repeated for all grid points. 

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

  print *,myproc,  'totalresiduei = ',totalresiduei &
               ,'  totalresiduee = ',totalresiduee 

!#ifndef subtract_average_rho
!! force global charge neutrality
!  np_i=min(np_i,np_e); np_e=min(np_i,np_e)
!print *,'np_i = ',np_i
!print *,'np_e = ',np_e
!#endif

  call redistribute(rvi,np_i)
  call redistribute(rve,np_e)

  call output(0,.true.)
End Subroutine init_waveH

