!******************************************************************************
!             3D particle code: Taylor Green Vortex.
!                                Tulasi Parashar - 2014/04/23
!
!******************************************************************************

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

Subroutine init_TGV()
  Use pe_env
  Use partfield
  Implicit None

  Real(kind=drk) gx,gy,gz,yy,b,rand_num,nsup,j,xx,zz
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez),&
                          pi2=6.2831853, pi4=2.*pi2, kval=pi2/lx
  Integer i,x,y,z,pe
  Integer :: numpart,numparti,numparte,numpartgrid
  Real(kind=drk) :: xp,yp,zp,xm,ym,zm
  Real(kind=drk) :: ni,ne,vix,vex&
                   ,viy,viz,vey,vez&
                   ,bxi,byi,bzi&
                   ,inter_vx, inter_vy, inter_vz&
                   ,inter_jx, inter_jy, inter_jz
  Real(kind=drk) :: sumni,sumne,globalsumni,globalsumne
  Real(kind=drk), dimension(nx,ny,nz) :: curl_x, curl_y, curl_z

! print *
!!!real(kind=drk) :: valtemp(nx)



! Runs from 0 to lx basically over the whole simulation.
  gx(x)=((x-0.5)+my_pex*nx)*dx
  gy(y)=((y-0.5)+my_pey*ny)*dy
  gz(z)=((z-0.5)+my_pez*nz)*dz

  ni(xx)     = 1.0
  ne(xx)     = 1.0

  bxi(xx,yy,zz)    =  bpert*cos(pi2*xx/lx)*sin(pi2*yy/ly)*sin(pi2*zz/lz)
  byi(xx,yy,zz)    =  bpert*sin(pi2*xx/lx)*cos(pi2*yy/ly)*sin(pi2*zz/lz)
  bzi(xx,yy,zz)    = -bpert*sin(pi2*xx/lx)*sin(pi2*yy/ly)*cos(pi2*zz/lz)*2.0

  vix(xx,yy,zz)    =  vpert*sin(pi2*xx/lx)*cos(pi2*yy/ly)*cos(pi2*zz/lz) 
  viy(xx,yy,zz)    = -vpert*cos(pi2*xx/lx)*sin(pi2*yy/ly)*cos(pi2*zz/lz)
  viz(xx,yy,zz)    =  0.

  vex(xx,yy,zz)    = vix(xx,yy,zz)
  vey(xx,yy,zz)    = viy(xx,yy,zz)
  vez(xx,yy,zz)    = viz(xx,yy,zz)

  nsup = 1.

! Anything processor dependent must be done after this next line.
  call init_pe_env()
  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** init: compiled for different number of PEs *****'
    call exitallpes()
  endif

!!!  print *,'myproc = ',myproc,'  vix = ',valtemp

  do z=1,nz ; do y=1,ny ; do x=1,nx
    sumni = sumni + ni(gx(x))
    sumne = sumne + ne(gx(x))
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
#ifdef relativistic
    write(6,*) '     relativistic initialization'
#endif
    write(6,*) '********** parameters **********'
  endif

! set electric and magnetic field (no kinetic equilibrium with guide field bz)
  do z=1,nz; do y=1,ny; do x=1,nx
    b1x(x,y,z)=bxi(gx(x),gy(y),gz(z))+b0x; 
    b1y(x,y,z)=byi(gx(x),gy(y),gz(z))+b0y;
    b1z(x,y,z)=bzi(gx(x),gy(y),gz(z))+b0z
    e1x(x,y,z)=0.; e1y(x,y,z)=0.; e1z(x,y,z)=0.
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
        vix(real(rvi(1,np_i),drk),real(rvi(2,np_i),drk),real(rvi(3,np_i),drk)),&
        viy(real(rvi(1,np_i),drk),real(rvi(2,np_i),drk),real(rvi(3,np_i),drk)),&
        viz(real(rvi(1,np_i),drk),real(rvi(2,np_i),drk),real(rvi(3,np_i),drk)),&
        )
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
  end do
end do ; end do ; end do
   do i = 1,np_e
    xx=int(rve(1,i)/dx-my_pex*nx+0.5)
    yy=int(rve(2,i)/dy-my_pey*ny+0.5)
    zz=int(rve(3,i)/dz-my_pez*nz+0.5)
    xm=rve(1,i)/dx-my_pex*nx+0.5-xx; xp=1.-xm
    ym=rve(2,i)/dy-my_pey*ny+0.5-yy; yp=1.-ym
    zm=rve(3,i)/dz-my_pez*nz+0.5-zz; zp=1.-zm
! interpolate the currents and then calculate the electron viscosity
! ve = vi - j/n (with q= +- 1)
    inter_jx = (curl_x(xx,yy,zz)*xp*yp*zp + &
      curl_x(xx,yy,zz+1)*xp*yp*zm + &
      curl_x(xx+1,yy,zz)*xm*yp*zp + &
      curl_x(xx+1,yy,zz+1)*xm*yp*zm + &
      curl_x(xx,yy+1,zz)*xp*ym*zp + &
      curl_x(xx,yy+1,zz+1)*xp*ym*zm + &
      curl_x(xx+1,yy+1,zz)*xm*ym*zp + &
      curl_x(xx+1,yy+1,zz+1)*xm*ym*zm)
    inter_jy = (curl_y(xx,yy,zz)*xp*yp*zp + &
      curl_y(xx,yy,zz+1)*xp*yp*zm + &
      curl_y(xx+1,yy,zz)*xm*yp*zp + &
      curl_y(xx+1,yy,zz+1)*xm*yp*zm + &
      curl_y(xx,yy+1,zz)*xp*ym*zp + &
      curl_y(xx,yy+1,zz+1)*xp*ym*zm + &
      curl_y(xx+1,yy+1,zz)*xm*ym*zp + &
      curl_y(xx+1,yy+1,zz+1)*xm*ym*zm)
    inter_jz = (curl_z(xx,yy,zz)*xp*yp*zp + &
      curl_z(xx,yy,zz+1)*xp*yp*zm + &
      curl_z(xx+1,yy,zz)*xm*yp*zp + &
      curl_z(xx+1,yy,zz+1)*xm*yp*zm + &
      curl_z(xx,yy+1,zz)*xp*ym*zp + &
      curl_z(xx,yy+1,zz+1)*xp*ym*zm + &
      curl_z(xx+1,yy+1,zz)*xm*ym*zp + &
      curl_z(xx+1,yy+1,zz+1)*xm*ym*zm)

      inter_vx = vix(real(rve(1,np_e),drk),real(rve(2,np_e),drk),real(rve(3,np_e),drk))-inter_jx/n_0 
      inter_vx = viy(real(rve(1,np_e),drk),real(rve(2,np_e),drk),real(rve(3,np_e),drk))-inter_jy/n_0 
      inter_vx = viz(real(rve(1,np_e),drk),real(rve(2,np_e),drk),real(rve(3,np_e),drk))-inter_jz/n_0 

     call veloadd(rve(4,i),rve(5,i),rve(6,i), &
       inter_vx, inter_vy, inter_vz)
   enddo

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
End Subroutine init_TGV

