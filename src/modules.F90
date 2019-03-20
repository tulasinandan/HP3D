!******************************************************************************
!                   3D particle code: module definition
!                           Andreas Zeiler, 1999
!
!                        LATEST CHANGE: May 9, 2001
!
!******************************************************************************

#include "param"

#ifndef hybrid
#ifndef twofluid
#define fullparticle
#endif
#endif

#ifndef RANDSEEDADD
#define RANDSEEDADD 54321
#endif

!------------------------------------------------------------------------------
!           standard PE related information, global sum/min/max
!------------------------------------------------------------------------------

Module pe_env
  Implicit	None

  include "mpif.h"

  Integer, Save :: myproc, nprocs
  Integer, Save :: mpi_drk, mpi_prk

  Integer,Save	:: my_pex, my_pey, my_pez
  Integer, Dimension (-1:+1,-1:+1,-1:+1), Save :: nb_pe
  Integer, Parameter :: n_pes=pex*pey*pez

  Integer, Dimension(8), Save :: alt_nb_2d
  Integer, Dimension(26), Save :: alt_nb_3d

  Logical, Save :: exchange_z_boundary

! my_pex, my_pey, my_pez: index of my_pe in x,y,z direction
! nb_pe(x,y,z): PE number of next neighbour (periodic boundaries)

Contains
  Subroutine init_pe_env
    use ran_state
    Implicit None
    Integer :: s_pex, s_pey, s_pez, x, y, z
    Integer :: mpi_err
    Integer :: seedarray, k

    s_pex(x) = Mod(x+pex,pex)
    s_pey(y) = Mod(y+pey,pey)
    s_pez(z) = Mod(z+pez,pez)

!   Set up MPI

    Call MPI_Init(mpi_err)
    if(mpi_err /= 0) then
      print *,'Error calling MPI_Init'
      stop
    endif
    call MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_ARE_FATAL,mpi_err)
    call MPI_Comm_size(MPI_COMM_WORLD,nprocs,mpi_err)
    call MPI_Comm_rank(MPI_COMM_WORLD,myproc,mpi_err)

    if (drk==4) then
      mpi_drk = MPI_REAL4
    else if(drk==8) then
      mpi_drk = MPI_REAL8
    else
      if(myproc==0) write(6,*) 'Default real kind = ',drk,' unsupported'
      stop
    endif

    if (prk==4) then
      mpi_prk = MPI_REAL4
    else if(prk==8) then
      mpi_prk = MPI_REAL8
    else
      if(myproc==0) write(6,*) 'Particle real kind = ',prk,' unsupported'
      stop
    endif

    my_pex = Mod(myproc,pex)
    my_pey = Mod(myproc,pex*pey)/pex
    my_pez = myproc/(pex*pey)
    do x=-1,1; do y=-1,1; do z=-1,1
      nb_pe(x,y,z) = &
             s_pex(my_pex+x)+pex*s_pey(my_pey+y)+(pex*pey)*s_pez(my_pez+z)
    enddo; enddo; enddo

!   Nearest neighbors in 1-D, not including the processor itself
    alt_nb_2d = (/ reshape(nb_pe(:,-1,0),(/3/)),nb_pe(-1,0,0),nb_pe(1,0,0),&
	reshape(nb_pe(:,1,0),(/3/)) /)
    alt_nb_3d = (/ reshape(nb_pe(:,:,-1),(/9/)),alt_nb_2d(:),&
	reshape(nb_pe(:,:,1),(/9/)) /)

    exchange_z_boundary = .true.

!   Set the random seed so that every PE will get a different random number 
!       sequence
!   RJ: I dont know if just setting the seed value to myproc will do that also,
!       the following formula is completely arbitrary

    seedarray = myproc*12345+myproc*myproc+RANDSEEDADD
    call ran_seed(sequence=seedarray)
!   seedarray(:) = myproc*12345+myproc*myproc+54321
!   call random_seed(size=k)
!   call random_seed(put=seedarray(1:k))

    return
  End Subroutine init_pe_env

  Function totalsum(sum)
    Implicit None
    Real(kind=drk) :: sum, sumtotal, totalsum
    Integer :: mpi_err
    call MPI_Allreduce(sum,sumtotal,1,mpi_drk,MPI_SUM,MPI_COMM_WORLD,mpi_err)
    totalsum=sumtotal
  End Function totalsum

  Function totalsumint(sum)
    Implicit None
    Integer :: sum, sumtotal, totalsumint
    Integer :: mpi_err
    call MPI_Allreduce(sum,sumtotal,1,mpi_integer,MPI_SUM,MPI_COMM_WORLD,mpi_err)
    totalsumint=sumtotal
  End Function totalsumint


  Function totalmax(data)
    Implicit None
    Real(kind=drk) :: data, maxtotal, totalmax
    Integer :: mpi_err
    call MPI_Allreduce(data,maxtotal,1,mpi_drk,MPI_MAX,MPI_COMM_WORLD,mpi_err)
    totalmax=maxtotal
  End Function totalmax
 
  Function reducesum(sum)
   ! MS; Reduction on processor 0
   Implicit None
   Real(kind=drk) :: sum, sumreduce, reducesum
   Integer :: mpi_err
   call MPI_Reduce(sum,sumreduce,1,mpi_drk,MPI_SUM,0,MPI_COMM_WORLD,mpi_err)
   reducesum=sumreduce
  End Function reducesum

  Function reducesum_logical(sum)
   ! MS; Reduction on processor 0
   Implicit None
   Logical :: sum, sumreduce, reducesum_logical
   Integer :: mpi_err
   call MPI_Reduce(sum,sumreduce,1,mpi_logical,MPI_LOR,0,MPI_COMM_WORLD,mpi_err)
   reducesum_logical=sumreduce
  End Function reducesum_logical 

  Function totalmin(data)
    Implicit None
    Real(kind=drk) :: data, mintotal, totalmin
    Integer :: mpi_err
    call MPI_Allreduce(data,mintotal,1,mpi_drk,MPI_MIN,MPI_COMM_WORLD,mpi_err)
    totalmin=mintotal
  End Function totalmin

  Subroutine exitallpes()
    Implicit None
    Integer mpi_err
    Call MPI_Abort(MPI_COMM_WORLD,1,mpi_err)
  End Subroutine exitallpes
 
End Module pe_env

!------------------------------------------------------------------------------
!                     global data: fields and particles
!------------------------------------------------------------------------------

Module partfield
  Implicit None

#ifdef fullparticle
! run as full particle code
  Real(kind=drk), Save, Dimension(0:nx+1,0:ny+1,0:nz+1) :: &
    e1x,e1y,e1z,b1x,b1y,b1z,e2x,e2y,e2z,b2x,b2y,b2z,esx,esy,esz,jx,jy,jz,rho &
   ,pxx,pyy,pzz,pxz,pyz,pxy
!     electric and magnetic field ([eb]..), smoothed electric field (es.)
!     at two different times (.[12].), three components (..[xyz]), 
!     j., rho current and charge density
!     pxx,pyy,pzz,pxz,pyz,pxy: Pressure tensors
  Real(kind=prk), Save, Dimension(6,maxparticles) :: rvi, rve
!     electron (..e) and ion (..i) particle coordinates (data fields 1-3)
!     and velocities (data fields 4-6)
  Integer, Save :: np_e, np_i
!     number of particles on PE
  Real(kind=drk), Save :: t, n_avg
! global time, number of particles per grid point (one species) according to
! normalization (n_avg=n0*dx*dy*dz)
#ifdef outputdist
  Real(kind=drk), Save :: disti(nx,nvxi), diste(nx,nvxe)
! Distribution functions of ions and electrons.
#endif
#endif
#ifdef hybrid
! run as hybrid code (particle ions, fluid electrons)
  Real(kind=drk), Save, Dimension(0:nx+1,0:ny+1,0:nz+1) :: &
    b1x,b1y,b1z,b2x,b2y,b2z,e1x,e1y,e1z,esx,esy,esz,pe1,pe2,jx,jy,jz,rho &
   ,pxx,pyy,pzz,pxz,pyz,pxy,jex,jey,jez,jtotx,jtoty,jtotz
! added jex jey jez on 11/13/2006 --- Tulasi Nandan
!     magnetic field (b..), and electron pressure (pe.)
!     at two different times, three components (..[xyz]), 
!     e1x,e1y,e1z is used as array for Eprime and E" (temporary)
!     esx,esy,esz is used for the smoothed electric field and as
!     temporary array for Bprime
!     jx,jy,jz is the ion current, and rho the ion density
!     (variable names had to be kept for compatibility with particle code)
  Real(kind=prk), Save, Dimension(6,maxparticles) :: rvi
! Real(kind=prk), Save, Dimension(6,maxtestparticles) :: rv1
!     ion particle coordinates (data fields 1-3)
!      and velocities (data fields 4-6)
  Integer, Save :: np_i
! Integer, Save :: np_1
!     number of particles on PE
  Real(kind=drk), Save :: t, n_avg
! global time, number of particles per grid point (one species) according to
! normalization (n_avg=n0*dx*dy*dz)
#endif

#ifdef twofluid
! run as two-fluid code
  Real(kind=drk), Save, Dimension(0:nx+1,0:ny+1,0:nz+1) :: &
    b1x,b1y,b1z,b2x,b2y,b2z,pe1,pe2,pi1,pi2,n1,n2,j1x,j1y,j1z,j2x,j2y,j2z
!     magnetic field (b..), pressure (p..), density (n.), and ion current (j..)
!     at two different times, three components (..[xyz]), 
!     e1x,e1y,e1z is used as array for Eprime and E" (temporary)
  Real(kind=drk), Save :: t,n_avg
! global time, dummy
#endif

#ifdef DRIVE_B
  Real(kind=drk), Save, Dimension(1:nx,1:ny,1:nz) :: dr_bx, dr_by, dr_bz
#endif
#ifdef DRIVE_V
  Real(kind=drk), Save, Dimension(1:nx,1:ny,1:nz) :: dr_vx, dr_vy, dr_vz
#endif

End Module partfield

#ifdef EXPANSION
!------------------------------------------------------------------------------
!                           Expansion parameters
!                                                 - TNP 12/12/12
!------------------------------------------------------------------------------
Module expansion_params
    
   Implicit None
   Real(kind=drk) :: ww, dwow, aa, daoa, usw
#ifdef ACCELERATED_EXPANSION
   Real(kind=drk) :: ti,ai,daoai,wi,dwowi,usi
   Real(kind=drk) :: tf,af,daoaf,wf,dwowf,usf
   Real(kind=drk) :: dtwd
#endif

End Module expansion_params
#endif

!------------------------------------------------------------------------------
!                     run time statistics of subroutines
!------------------------------------------------------------------------------

Module timer

  Use pe_env

  Implicit none

  Real(kind=8) :: time_start, time_stop
  Real(kind=drk) :: t_stepfield, t_stepx, t_calc_rho, t_smooth_rho, &
                    t_divergence, t_stepv, t_calc_j, t_smooth_j, t_energy, &
                    t_movie, t_partsort, t_stephybrid, t_calc_e, t_step2fluid,&
		    t_redis,t_calc_p,t_calc_je
  Integer :: n_stepfield, n_stepx, n_calc_rho, n_smooth_rho, n_divergence, &
             n_stepv, n_calc_j, n_smooth_j, n_energy, n_movie, n_partsort, &
             n_stephybrid, n_calc_e, n_step2fluid,n_redis,n_calc_p,n_calc_je
!------------------------------------------
! t_calc_je, n_calc_je added by Tulasi
! on 11/13/2006
!------------------------------------------

Contains
  Subroutine init_timer()
    Implicit None
    t_stepfield=0.; n_stepfield=0
    t_stepx=0.; n_stepx=0
    t_calc_rho=0.; n_calc_rho=0
    t_smooth_rho=0.; n_smooth_rho=0
    t_divergence=0.; n_divergence=0
    t_stepv=0.; n_stepv=0
    t_calc_j=0.; n_calc_j=0
    t_smooth_j=0.; n_smooth_j=0
    t_energy=0.; n_energy=0
    t_movie=0.; n_movie=0
    t_partsort=0.; n_partsort=0
    t_stephybrid=0.; n_stephybrid=0
    t_calc_e=0.; n_calc_e=0
    t_step2fluid=0.; n_step2fluid=0
    t_redis = 0.; n_redis = 0
    t_calc_p = 0.; n_calc_p = 0
! t_calc_je, n_calc_je added by Tulasi
! on 11/13/2006
    t_calc_je = 0.; n_calc_je=0
  End Subroutine init_timer

  Subroutine start_timer()
    Implicit None
    time_start = MPI_Wtime()
  End Subroutine start_timer

  Subroutine stop_timer(t,n)
    Implicit None
    Real(kind=drk) :: t
    Integer :: n
    time_stop = MPI_Wtime()
    t=t+(time_stop-time_start)
    n=n+1
  End Subroutine stop_timer

  Subroutine read_timer()
    Implicit None
    Real(kind=drk) :: total_time
	total_time = t_stepfield+t_stepx+t_calc_rho+t_smooth_rho+ &
                     t_divergence+t_stepv+t_calc_j+t_smooth_j+t_energy+&
                     t_movie+t_partsort+t_stephybrid+t_calc_e+t_step2fluid+&
		     t_redis+t_calc_p
    if (myproc==0) then
	  write(6,*) '************* performance **************'
      write(6,*) '**** average run time of subroutines ***'
      if (n_stepfield .ne. 0) &
        write(6,"(A,F9.5,A)") '  stepfield:  ',t_stepfield/n_stepfield,' s'
      if (n_stephybrid .ne. 0) &
        write(6,"(A,F9.5,A)") '  stephybrid: ',t_stephybrid/n_stephybrid,' s'
      if (n_step2fluid .ne. 0) &
        write(6,"(A,F9.5,A)") '  step2fluid: ',t_step2fluid/n_step2fluid,' s'
      if (n_stepx .ne. 0) &
        write(6,"(A,F9.5,A)") '  stepx:      ',t_stepx/n_stepx,' s'
      if (n_calc_rho .ne. 0) &
        write(6,"(A,F9.5,A)") '  calc_rho:   ',t_calc_rho/n_calc_rho,' s'
      if (n_smooth_rho .ne. 0) &
        write(6,"(A,F9.5,A)") '  smooth_rho: ',t_smooth_rho/n_smooth_rho,' s'
      if (n_divergence .ne. 0) &
        write(6,"(A,F9.5,A)") '  divergence: ',t_divergence/n_divergence,' s'
      if (n_stepv .ne. 0) &
        write(6,"(A,F9.5,A)") '  stepv:      ',t_stepv/n_stepv,' s'
      if (n_calc_j .ne. 0) &
        write(6,"(A,F9.5,A)") '  calc_j:     ',t_calc_j/n_calc_j,' s'
      if (n_smooth_j .ne. 0) &
        write(6,"(A,F9.5,A)") '  smooth_j:   ',t_smooth_j/n_smooth_j,' s'
      if (n_calc_p .ne. 0) &
        write(6,"(A,F9.5,A)") '  calc_p:   ',t_calc_p/n_calc_p,' s'
      if (n_calc_e .ne. 0) &
        write(6,"(A,F9.5,A)") '  calc_e:     ',t_calc_e/n_calc_e,' s'
      if (n_energy .ne. 0) &
        write(6,"(A,F9.5,A)") '  energy:     ',t_energy/n_energy,' s'
      if (n_movie .ne. 0) &
        write(6,"(A,F9.5,A)") '  movie:      ',t_movie/n_movie,' s'
      if (n_partsort .ne. 0) &
        write(6,"(A,F9.5,A)") '  partsort:   ',t_partsort/n_partsort,' s'
      if (n_redis .ne. 0) &
        write(6,"(A,F9.5,A)") '  redistrib.: ',t_redis/n_redis,' s'

      if (total_time .ne. 0.) then
        write(6,*) '**** time spent in subroutines ***'
	    write(6,"(A,F9.1,A)") '  diagnosed time: ', total_time,' s'
        write(6,"(A,F11.2,A)") &
                  '  stepfield:  ',t_stepfield/total_time*100.,' %'
        write(6,"(A,F11.2,A)") '  stepx:      ',t_stepx/total_time*100.,' %'
        write(6,"(A,F11.2,A)") '  calc_rho:   ',t_calc_rho/total_time*100.,' %'
        write(6,"(A,F11.2,A)") & 
                  '  smooth_rho: ',t_smooth_rho/total_time*100.,' %'
        write(6,"(A,F11.2,A)") &
                  '  divergence: ',t_divergence/total_time*100.,' %'
        write(6,"(A,F11.2,A)") '  stepv:      ',t_stepv/total_time*100.,' %'
        write(6,"(A,F11.2,A)") '  calc_j:     ',t_calc_j/total_time*100.,' %'
        write(6,"(A,F11.2,A)") '  calc_p:     ',t_calc_p/total_time*100.,' %'
        write(6,"(A,F11.2,A)") '  smooth_j:   ',t_smooth_j/total_time*100.,' %'
        write(6,"(A,F11.2,A)") '  calc_e:     ',t_calc_e/total_time*100.,' %'
        write(6,"(A,F11.2,A)") '  energy:     ',t_energy/total_time*100.,' %'
        write(6,"(A,F11.2,A)") '  movie:      ',t_movie/total_time*100.,' %'
        write(6,"(A,F11.2,A)") '  partsort:   ',t_partsort/total_time*100.,' %'
        write(6,"(A,F11.2,A)") '  redistrib.: ',t_redis/total_time*100.,' %'
        write(6,"(A,F11.2,A)") &
                  '  stephybrid: ',t_stephybrid/total_time*100.,' %'
        write(6,"(A,F11.2,A)") &
                  '  step2fluid: ',t_step2fluid/total_time*100.,' %'
	    write(6,*) '************* performance **************'
      endif
    endif
  End Subroutine read_timer
End Module timer

!#ifndef CRAY
!Subroutine flush(lu)
!  Implicit None
!  Integer :: lu
!
!! Dummy for Non-Cray machines
!
!End Subroutine flush
!#endif
