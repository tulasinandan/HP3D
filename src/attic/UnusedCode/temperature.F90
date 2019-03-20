!******************************************************************************
!                   3D particle code: velocity distribution
!                           Andreas Zeiler, 1999
!
!                        LATEST CHANGE: May 3, 2001
!
!******************************************************************************

#include "param"

Module channels
  Implicit  None
  Integer, Dimension(nchannels) :: unit
  Integer :: ngridx,ngridy,ngridz
  Integer :: count
  Character*16 :: filename
End Module channels

!------------------------------------------------------------------------------

Program temperature
  Use pe_env
  Use channels
  Implicit None

  Character*16 :: argbuf
  Integer :: mpi_err

  call init_pe_env()

! read filenumber, gridsize
  call our_getarg(1,argbuf)
  if (argbuf == ' ') then
    write(6,*) '***** temperature: missing number of input file *****'
    stop
  endif
  read(argbuf,*) count
  call our_getarg(2,argbuf); if (argbuf.ne.' ') read(argbuf,*) ngridx
  call our_getarg(3,argbuf); if (argbuf.ne.' ') read(argbuf,*) ngridy
  call our_getarg(4,argbuf); if (argbuf.ne.' ') read(argbuf,*) ngridz

  if (myproc == 0) then
    write(filename,"(a,'-xxx.',i3.3)") head, count
    write(6,*) 'input from ', filename
    write(6,*) '    ','ngridx = ',ngridx
    write(6,*) '    ','ngridy = ',ngridy
    write(6,*) '    ','ngridz = ',ngridz
  endif

  call main_prog()

  call MPI_Finalize(mpi_err)
End Program temperature

!------------------------------------------------------------------------------

Subroutine main_prog()
  Use pe_env
  Use channels
  Implicit None

  Real(kind=prk), Dimension(ngridx,ngridy,ngridz,3) :: sumvavg_ion,sumvavg_ele,vavg_ion,vavg_ele
  Real(kind=prk), Dimension(ngridx,ngridy,ngridz,6) :: sumv2_ion,sumv2_ele,v2_ion,v2_ele
  Integer, Dimension(ngridx,ngridy,ngridz) ::  sumnions,sumneles,nions,neles
  Integer :: n1,n2,n3,n4,n5,ch,i
  Integer :: mpi_err,loopcount

  Real(kind=drk) :: t,n_avg

! For large values of ngridx/y/z this routine can cause a segmentation
! fault.  The problem, I think, is a memory overflow because of the
! many large arrays.  Perhaps smarter memory management is in order.

  sumnions = 0.; sumneles = 0.; sumv2_ion = 0.; sumv2_ele = 0.
  sumvavg_ion = 0.; sumvavg_ele = 0.

  do loopcount = 1,ceiling(nchannels*1./nprocs)
    ch = (myproc+1) + (loopcount-1)*nprocs
    if (ch .gt. nchannels) exit 
    unit(ch)=30+ch
    write(filename,"(a,'-',i3.3,'.',i3.3)") head, ch, count
    open(unit=unit(ch),file=filename,action='read',form='unformatted')
    read(unit(ch)) t,n_avg,n1,n2,n3,n4,n5
    if (n1 .ne. nx*pex .or. n2 .ne. ny*pey .or. n3 .ne. nz*pez &
                     .or. n4 .ne. bufsize .or. n5 .ne. nchannels) then
      write(6,*) '***** temperature: input-file incompatible *****'
 write(6,*) n1,nx,pex,n2,ny,pey,n3,nz,pez,n4,bufsize,n5,nchannels
      call MPI_Finalize(mpi_err)
      stop
    endif

! skip field data
    if (ch == 1) then
      call skipfield()
    endif

! read and diagnose ions
    call diagparticles(nions,v2_ion,vavg_ion,ch)
    sumnions = sumnions+nions; sumv2_ion = sumv2_ion+v2_ion
    sumvavg_ion = sumvavg_ion+vavg_ion
#ifndef hybrid
! read and diagnose electrons
    call diagparticles(neles,v2_ele,vavg_ele,ch)
    sumneles = sumneles+neles; sumv2_ele = sumv2_ele+v2_ele
    sumvavg_ele = sumvavg_ele+vavg_ele
#endif

    close(unit(ch))
  enddo

! Combine data. Note that I get tricky here by reusing vix, viy, vex, etc.
! Here they become the _result_ of summing sumvix over all processors.
  call MPI_Allreduce(sumnions,nions,ngridx*ngridy*ngridz,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(sumv2_ion,v2_ion,ngridx*ngridy*ngridz*6,mpi_prk,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(sumvavg_ion,vavg_ion,ngridx*ngridy*ngridz*3,mpi_prk,MPI_SUM,MPI_COMM_WORLD,mpi_err)
#ifndef hybrid
  call MPI_Allreduce(sumneles,neles,ngridx*ngridy*ngridz,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(sumv2_ele,v2_ele,ngridx*ngridy*ngridz*6,mpi_prk,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(sumvavg_ele,vavg_ele,ngridx*ngridy*ngridz*3,mpi_prk,MPI_SUM,MPI_COMM_WORLD,mpi_err)
#endif

  do i = 1,3
    vavg_ion(:,:,:,i) = vavg_ion(:,:,:,i)/nions
    v2_ion(:,:,:,i) = v2_ion(:,:,:,i)/nions
    vavg_ele(:,:,:,i) = vavg_ele(:,:,:,i)/neles
    v2_ele(:,:,:,i) = v2_ele(:,:,:,i)/neles
  enddo

  do i = 4,6
    v2_ion(:,:,:,i) = v2_ion(:,:,:,i)/nions
    v2_ele(:,:,:,i) = v2_ele(:,:,:,i)/neles
  enddo

! Temperatures have three components.  Multiply by electron mass here, but
! do not sum over components; leave that for IDL. Reuse sumv2_ion/ele to
! mean ion_t/ele_t.
  sumv2_ion(:,:,:,1:3) = v2_ion(:,:,:,1:3) - vavg_ion**2
  sumv2_ion(:,:,:,4) = v2_ion(:,:,:,4) - vavg_ion(:,:,:,1)*vavg_ion(:,:,:,2)
  sumv2_ion(:,:,:,5) = v2_ion(:,:,:,5) - vavg_ion(:,:,:,2)*vavg_ion(:,:,:,3)
  sumv2_ion(:,:,:,6) = v2_ion(:,:,:,6) - vavg_ion(:,:,:,1)*vavg_ion(:,:,:,3)

  sumv2_ele(:,:,:,1:3) = v2_ele(:,:,:,1:3) - vavg_ele**2
  sumv2_ele(:,:,:,4) = v2_ele(:,:,:,4) - vavg_ele(:,:,:,1)*vavg_ele(:,:,:,2)
  sumv2_ele(:,:,:,5) = v2_ele(:,:,:,5) - vavg_ele(:,:,:,2)*vavg_ele(:,:,:,3)
  sumv2_ele(:,:,:,6) = v2_ele(:,:,:,6) - vavg_ele(:,:,:,1)*vavg_ele(:,:,:,3)
  sumv2_ele = m_e*sumv2_ele

  if (myproc == 0) then
    write(6,*)
    write(6,"(A)") '** Ion temperature **'
    write(6,*) sumv2_ion
#ifndef hybrid
    write(6,*)
    write(6,"(A)") '** Electron temperature **'
    write(6,*) sumv2_ele
#endif
  endif

  return
End Subroutine main_prog


!------------------------------------------------------------------------------

Subroutine skipfield()
  Use channels
  Implicit None
  Integer x,y,z,px,py,pz,ch,fields
  Real(kind=drk), Dimension(1:nx*pex) :: tmp

! read E and B fields from original files and discard them (no use here)
#ifndef hybrid
  do fields=1,6
#else
  do fields=1,4
#endif
    do pz=0,pez-1; do z=1,nz; do py=0,pey-1; do y=1,ny
      ch=mod(pz*nz+z-1,nchannels)+1 ! New I/O Scheme
      read(unit(ch)) tmp
    enddo; enddo; enddo; enddo
  enddo

  return
End Subroutine skipfield

!------------------------------------------------------------------------------

Subroutine diagparticles(npartarr,v2_tot_arr,vavg_arr,inchannel)
  Use channels
  Implicit None
  Real(kind=prk), Dimension(ngridx,ngridy,ngridz,3) :: vavg_arr
  Real(kind=prk), Dimension(ngridx,ngridy,ngridz,6) :: v2_tot_arr
  Integer, Dimension(ngridx,ngridy,ngridz) :: npartarr
  Real(kind=prk), Dimension(6,maxparticles), Save :: rv
  Real(kind=prk), Dimension(6,bufsize) :: buffer
  Integer :: np,pos,n_pes_in,ch,i,j
  Integer :: xcell,ycell,zcell
  Integer :: inchannel,fudge

! read number of original PEs
  read(unit(inchannel)) n_pes_in

npartarr = 0
v2_tot_arr = 0.
vavg_arr = 0.

! Need to do something a little fancy to make sure we read the correct
! amount of data from each dump file.
  fudge = 0
  if (inchannel .le. mod(n_pes_in,nchannels)) fudge = 1

  do j = 0, n_pes_in/nchannels - 1 + fudge
! read particle buffer of all PEs in given dump file
    read(unit(inchannel)) np
    pos=0
    do while (pos<np)
      read(unit(inchannel)) buffer
      rv(:,pos+1:pos+bufsize) = buffer(:,:)
      pos=pos+bufsize
        if (pos+bufsize .gt. maxparticles) then
          write(6,*) '***** vdist_par: particle buffer overflow *****'
          stop
        endif
    enddo
! analyse particle buffer
    do i=1,np
! find x,y,z cells 
xcell = int(rv(1,i)/(lx/ngridx))+1
ycell = int(rv(2,i)/(ly/ngridy))+1
zcell = int(rv(3,i)/(lz/ngridz))+1

npartarr(xcell,ycell,zcell) = npartarr(xcell,ycell,zcell) + 1

! should do this as array total
v2_tot_arr(xcell,ycell,zcell,1:3) = v2_tot_arr(xcell,ycell,zcell,1:3) + rv(4:6,i)**2
v2_tot_arr(xcell,ycell,zcell,4) = v2_tot_arr(xcell,ycell,zcell,4) + rv(4,i)*rv(5,i)
v2_tot_arr(xcell,ycell,zcell,5) = v2_tot_arr(xcell,ycell,zcell,5) + rv(5,i)*rv(6,i)
v2_tot_arr(xcell,ycell,zcell,6) = v2_tot_arr(xcell,ycell,zcell,6) + rv(4,i)*rv(6,i)
vavg_arr(xcell,ycell,zcell,1:3) = vavg_arr(xcell,ycell,zcell,1:3) + rv(4:6,i)

    enddo
  enddo
  return 
End Subroutine diagparticles


Subroutine our_getarg(num,buf)

  Implicit none
  Integer num
  Character*(*) buf
#ifdef CRAY
  Integer ilen, ierr
#endif

  buf = ' '

#ifdef CRAY
  call pxfgetarg(num,buf,ilen,ierr)
#endif
#ifdef AIX
  call getarg(num,buf)
#endif
#ifdef HITACHI
  call getarg(num+1,buf)
#endif

End
