!================================================================================
!
!
!
!   This is the initialization subroutine for the hybrid code to understand energy 
!   conservation.
!
!
!		Tulasi Nandan Parashar 01/02/2009
!
!
!
!================================================================================

#include "param"
#ifndef b0x
#define b0x 0.
#endif
#ifndef b0y
#define b0y 0.
#endif
#ifndef b0z
#define b0z 0.
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

Subroutine init_anydist()
  Use pe_env
  Use partfield
  Implicit None

  Real(kind=drk) gx,gy,gz,xx,yy,zz,b,n,rand_num,sech,ppert,npert
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez),&
                               pi2=6.2831853, pi4=2*pi2
  Integer i,x,y,z,pe,sign,nsup,distsize, vi, vj, vk
  Real(kind=drk) :: sum_ncell, ntot_p, vii, vjj, vkk
  Integer :: mpi_err
  Integer :: part_counter
  Integer, Parameter :: nbinsx=201, nbinsy=201, nbinsz=201
  Real(kind=drk), Parameter :: vxmin=-2.0, vxmax=2.0, vymin=-2.0, vymax=2.0,&
                               vzmin=-2.0, vzmax=2.0, vxd=vxmax-vxmin, vyd=vymax&
                               -vymin, vzd=vzmax-vzmin, dvx=vxd/nbinsx, dvy=vyd/&
                               nbinsy, dvz=vzd/nbinsz
  Real(kind=drk), Dimension(nbinsx,nbinsy,nbinsz) :: distfn
  gx(x)=((x-0.5)+my_pex*nx)*dx
  gy(y)=((y-0.5)+my_pey*ny)*dy
  gz(z)=((z-0.5)+my_pez*nz)*dz
  sech(yy)=1/cosh(yy)

  ppert(xx,yy)=0.

  npert(xx,yy)= ppert(xx,yy)/T_i
  nsup=int(1./(T_i+T_e)+5.)

  call init_pe_env()
  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** init-hybrid: compiled for different number of PEs *****'
    call exitallpes()
  endif

!write parameters to log-file
  if (myproc==0) then
    write(6,*) '********** initnone **********'
    write(6,*) '******* (hybrid version) *******'
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
    write(6,*) '     part./gridp. = ',ppg
    write(6,*) '     T_i = ',T_i
    write(6,*) '     T_e = ',T_e
#ifndef no_background_density
    write(6,*) '     n_0 = ',n_0
#else
    write(6,*) '     n_0 = 0'
#endif
    write(6,*) '     b0 = ',b0x, b0y, b0z
    write(6,*) '     w0 = ',w0
    write(6,*) '     psi0 = ',psi0
    write(6,*) '********** parameters **********'
  endif

!! Read in the distribution function
  inquire(iolength=distsize) distfn
  open(unit=931,file='distfn.dat',form='unformatted',access='direct',recl=distsize)
  read(931,rec=1) distfn
  close(931)

! set magnetic field
  do z=1,nz; do y=1,ny; do x=1,nx
    b1x(x,y,z)=b0x; b1y(x,y,z)=b0y; b1z(x,y,z)=b0z
  enddo; enddo; enddo 
!-------------------------------------------
! finding the total number of particles
! Ntot_proc=n_avg*Sum(n_cell)/(pex*pey)
!	 		Tulasi
  ntot_p=0
  do x=1,nx; do y=1,ny
    ntot_p = ntot_p + npert(gx(x),gy(y)) 
  enddo; enddo
  ntot_p = nz*ntot_p
  
  call MPI_Allreduce(ntot_p,sum_ncell,1,mpi_drk,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  n_avg=ppg/n_0 
!-------------------------------------------
  np_i=0
  part_counter=0
! load ions in current sheets
!
 inf_loop: do
 !old  do i=1,nx*ny*nz*ppg*nsup/n_0
    np_i=np_i+1
    if (np_i .gt. maxparticles) then
      write(6,*) '***** init: particle buffer overflow *****'
      call exitallpes()
 !old    endif
 ! This is to make sure that we get the number of particles
 ! we want. We did not believe the previous version.
    elseif (np_i .ge. ppg*nx*ny*nz+int(ntot_p*n_avg)) then
    exit inf_loop
    else
! Find the random location as well as the random velocity of the particle.
! Accept or reject particles with probability based on local density as
! well as the distribution function.
    call location(rvi(1,np_i),rvi(2,np_i),rvi(3,np_i))
    vii = rand_num()*nbinsx; rvi(4,np_i) = vxmin+vii*dvx; vi = floor(vii)+1
    vjj = rand_num()*nbinsy; rvi(5,np_i) = vymin+vjj*dvy; vj = floor(vjj)+1
    vkk = rand_num()*nbinsz; rvi(6,np_i) = vzmin+vkk*dvz; vk = floor(vkk)+1
!   rvi(4,np_i) = vxmin+rand_num()*vxd; vi=floor(rvi(4,np_i)/dvx)+1
!   rvi(5,np_i) = vymin+rand_num()*vyd; vj=floor(rvi(5,np_i)/dvy)+1
!   rvi(6,np_i) = vzmin+rand_num()*vzd; vk=floor(rvi(6,np_i)/dvz)+1
!   if ((vi .gt. nbinsx) .or. (vj .gt. nbinsy) .or. (vk .gt. nbinsz) ) then
!   write(6,*) 'Vel index out of range' 
!   exit inf_loop
!   endif

! If not in density and/or distfn range, discard it
    if ((rand_num()*nsup .ge. n_0+npert(real(rvi(1,np_i),drk),&
       real(rvi(2,np_i),drk))) .and. (rand_num() .gt. distfn(vi,vj,vk))) then

      np_i=np_i-1
    endif

    endif ! for the maxparticles if loop
 !old  enddo
    enddo inf_loop
!-------------------------------------------

   n_avg=ppg/n_0  
  call redistribute(rvi,np_i)

! calculate electron pressure
    call calc_rho(rvi,np_i,real(+1.,drk),.true.)
  do z=1,nz; do y=1,ny; do x=1,nx
#ifdef barotropic
    pe1(x,y,z)=T_e*rho(x,y,z)**gam                
#else
    pe1(x,y,z)=rho(x,y,z)*T_e
#endif
  enddo; enddo; enddo 
! write start-up file
  call output(0,.true.)
End Subroutine init_anydist
