!******************************************************************************
!                   3D particle code: input/output package
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

Module channels
  Implicit  None
  Integer, Dimension(nchannels) :: pe,unit
End Module channels

!------------------------------------------------------------------------------
!							   output
!------------------------------------------------------------------------------

Subroutine output(count,verbose)
  Use pe_env
  Use partfield
  Use channels
  Implicit None
  Character*16 :: filename
  Integer :: count,ch
  Logical :: verbose
  Real(kind=8) :: time1,time2
  Integer :: mpi_err


  call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
  time1 = MPI_Wtime()
  if (verbose .and. myproc == 0) then
    write(filename,"(a,'-xxx.',i3.3)") head, count
    write(6,*) 'time t=',t,'   output (kind=',prk,') to ', filename
  endif
! assign output channels to PEs and open files
  do ch=1,nchannels
    pe(ch)=mod(ch-1,n_pes)
    unit(ch)=30+ch
    if (myproc == pe(ch)) then
#ifdef USE_SIMPLE_IO  
      write(filename,"(a,'-',i4.4,'.',i3.3)") head, ch, count
#else
      write(filename,"(a,'-',i3.3,'.',i3.3)") head, ch, count
#endif
      open(unit=unit(ch),file=filename,action='write',form='unformatted')
      write(unit(ch)) t,n_avg,nx*pex,ny*pey,nz*pez,bufsize,nchannels
    endif
  enddo
! write data
#ifdef fullparticle
  call writefield(e1x); call writefield(e1y); call writefield(e1z)
  call writefield(b1x); call writefield(b1y); call writefield(b1z)
  call writeparticles(rvi,np_i); call writeparticles(rve,np_e)
#endif
#ifdef hybrid
  if (myproc == 0) write(*,*) 'Inside output, writing fields'
  call writefield(b1x); call writefield(b1y); call writefield(b1z)
  if (myproc == 0) write(*,*) 'Inside output, writing particles'
  call writefield(pe1); call writeparticles(rvi,np_i)
#endif
#ifdef twofluid
  call writefield(b1x); call writefield(b1y); call writefield(b1z)
  call writefield(pe1); call writefield(pi1); call writefield(n1)
  call writefield(j1x); call writefield(j1y); call writefield(j1z)
#endif
! close files
  do ch=1,nchannels
    if (myproc == pe(ch)) close(unit(ch))
  enddo
  call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
  time2 = MPI_Wtime()
  if (myproc == 0) write(6,*) 'output took ',(time2-time1),' seconds'
  return
End Subroutine output

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Subroutine writefield(f)
  Use pe_env
  Use channels
  Implicit None
  Integer x,y,z,px,py,pz,ch
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: f
  Real(kind=drk), Dimension(1:nx*pex) :: tmp
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status
  Integer :: mpi_err
  Integer npe

#ifdef USE_SIMPLE_IO

! 8/6/2008: M.Shay: Very simple I/O scheme

! New I/O scheme: slices of data are distributed across the output files

! Each processor outputs all of its data to a single file.

  ch=myproc+1
  write (unit(ch)) f(1:nx,1:ny,1:nz)

  call MPI_Barrier(MPI_COMM_WORLD,mpi_err)


#elif USE_OLD_IO

! Old I/O scheme: lines of data are distributed across the output files

  do pz=0,pez-1; do z=1,nz
    do py=0,pey-1; do y=1,ny
      ch=mod(py*ny+y-1,nchannels)+1
      if (myproc == pe(ch)) then
        do px=0,pex-1
          npe = px+pex*py+pex*pey*pz
          if (npe /= myproc) then
            call MPI_Recv(tmp(px*nx+1),nx,mpi_drk,npe,123,MPI_COMM_WORLD,mpi_status,mpi_err)
          else
            tmp(px*nx+1:px*nx+nx) = f(1:nx,y,z)
          endif
        enddo
        write (unit(ch)) tmp
      elseif (my_pez==pz .and. my_pey==py) then
!
!       Use MPI_Ssend here to avoid congestion on the receiving side
!
        call MPI_Ssend(f(1,y,z),nx,mpi_drk,pe(ch),123,MPI_COMM_WORLD,mpi_err)
      endif
    enddo; enddo
  enddo; enddo

#else

! New I/O scheme: slices of data are distributed across the output files

  Integer :: mpi_request
  Real(kind=drk), Dimension(0:nx+1,1:ny,0:pex-1) :: tmpr

! Every PE sends all of its slices to the receiving PE using MPI_Isend
! Important:
!   There may be a limit on outstanding MPI_Isend requests depending
!   on the MPI implementation. If nz is bigger than this limit, the
!   following code may fail !!!!
!   In this case the limit should be increased.
!   (On SGI/Cray systems by setting the environment variable MPI_REQUEST_MAX)

  do z=1,nz
    ch = mod(my_pez*nz+z-1,nchannels)+1
    call MPI_Isend(f(0,1,z),(nx+2)*ny,mpi_drk,pe(ch),123,MPI_COMM_WORLD,mpi_request,mpi_err)
    call MPI_Request_free(mpi_request,mpi_err)
  enddo

! The PEs responsible for output gather the slices and write them out

  do ch = 1,nchannels
    if (myproc == pe(ch)) then
      do z=ch,pez*nz,nchannels
        pz = (z-1)/nz
        do py=0,pey-1
          do px=0,pex-1
            npe = px+pex*py+pex*pey*pz
            call MPI_Recv(tmpr(0,1,px),(nx+2)*ny,mpi_drk,npe,123,MPI_COMM_WORLD,mpi_status,mpi_err)
          enddo
          do y=1,ny
            do px=0,pex-1
              tmp(px*nx+1:px*nx+nx) = tmpr(1:nx,y,px)
            enddo
            write (unit(ch)) tmp
          enddo
        enddo
      enddo
    endif
  enddo

  call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

#endif

  return
End Subroutine writefield

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifndef twofluid
Subroutine writeparticles(rv,np)
  Use pe_env
  Use channels
  Implicit None
  Real(kind=prk), Dimension(6,maxparticles) :: rv
  Real(kind=prk), Dimension(6,bufsize) :: buffer
  Integer :: npp
  Integer :: np,ppe,pos,i,j,ch
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status
  Integer :: mpi_err

#ifdef USE_SIMPLE_IO
! Very simple IO. Each processor just outputs what it has.

  ch=myproc+1
  write(unit(ch)) n_pes
  write(unit(ch)) np
  write(unit(ch)) rv(:,1:np)

#else
! Original input and output.

  do i=np+1,maxparticles; do j=1,6
    rv(j,i)=0.
  enddo; enddo
  do ch=1,nchannels
    if (myproc == pe(ch)) write(unit(ch)) n_pes
  enddo

! We use the fact that:  pe(ch)=mod(ch-1,n_pes)

  if(myproc < nchannels) then

    ! We have to write, first our own stuff, then the stuff
    ! from PEs myproc+nchannels, myproc+2*nchannels ...

    ch = myproc+1

    do ppe=myproc,n_pes-1,nchannels
      if(ppe /= myproc) then
        call MPI_Recv(npp,1,MPI_INTEGER,ppe,456,MPI_COMM_WORLD,mpi_status,mpi_err)
      else
        npp = np
      endif
      write(unit(ch)) npp
      pos=0
      do while (pos<npp)
        if(ppe /= myproc) then
          call MPI_Recv(buffer,6*bufsize,mpi_prk,ppe,456,MPI_COMM_WORLD,mpi_status,mpi_err)
        else
          buffer(:,:) = rv(:,pos+1:pos+bufsize)
        endif
        write(unit(ch)) buffer
        pos=pos+bufsize
      enddo
    enddo

  else

    ! We have to send our stuff to PE mod(myproc,nchannels)

    ! Use MPI_Ssend here to avoid congestion on the receiving side

    call MPI_Ssend(np,1,MPI_INTEGER,mod(myproc,nchannels),456,MPI_COMM_WORLD,mpi_err)
    pos=0
    do while (pos<np)
      call MPI_Ssend(rv(1,pos+1),6*bufsize,mpi_prk,mod(myproc,nchannels), &
                     456,MPI_COMM_WORLD,mpi_err)
      pos=pos+bufsize
    enddo

  endif
#endif

  return 
End Subroutine writeparticles
#endif

!------------------------------------------------------------------------------
!							   input
!------------------------------------------------------------------------------

Subroutine input(count,verbose)
  Use pe_env
  Use partfield
  Use channels
  Implicit None
  Character*16 :: filename
  Integer :: count,ppe,nn1,nn2,n3,n4,n5,ch
  Logical :: verbose
  Real(kind=8) :: time1,time2
  Integer mpi_err

  call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
  time1 = MPI_Wtime()
  if (verbose .and. myproc == 0) then 
    write(filename,"(a,'-xxx.',i3.3)") head, count
    write(6,*) 'output (kind=',prk,') from ', filename
  endif
! assign input channels to PEs and open files
  do ch=1,nchannels
    pe(ch)=mod(ch-1,n_pes)
    unit(ch)=30+ch
    if (myproc == pe(ch)) then
#ifdef USE_SIMPLE_IO  
      write(filename,"(a,'-',i4.4,'.',i3.3)") head, ch, count
#else
      write(filename,"(a,'-',i3.3,'.',i3.3)") head, ch, count
#endif
      open(unit=unit(ch),file=filename,action='read',form='unformatted')
      read(unit(ch)) t,n_avg,nn1,nn2,n3,n4,n5
      if (myproc == 0) then
        if (nn1 .ne. nx*pex .or. nn2 .ne. ny*pey .or. n3 .ne. nz*pez &
                           .or. n4 .ne. bufsize .or. n5 .ne. nchannels) then
          write(6,*) '***** input: input-file incompatible *****'
          call exitallpes()
        endif
      endif
    endif
  enddo
  call MPI_Bcast(t,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
  call MPI_Bcast(n_avg,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
! read data
#ifdef fullparticle
  call readfield(e1x); call readfield(e1y); call readfield(e1z)
  call readfield(b1x); call readfield(b1y); call readfield(b1z)
  call readparticles(rvi,np_i); call readparticles(rve,np_e)
#endif
#ifdef hybrid
  call readfield(b1x); call readfield(b1y); call readfield(b1z)
  call readfield(pe1); call readparticles(rvi,np_i)
#endif
#ifdef twofluid
  call readfield(b1x); call readfield(b1y); call readfield(b1z)
  call readfield(pe1); call readfield(pi1); call readfield(n1)
  call readfield(j1x); call readfield(j1y); call readfield(j1z)
#endif
! close files
  do ch=1,nchannels
    if (myproc == pe(ch)) close(unit(ch))
  enddo
  call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
  time2 = MPI_Wtime()
  if (myproc == 0) write(6,*) 'input took ',(time2-time1),' seconds'
  return
End Subroutine input

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Subroutine readfield(f)
  Use pe_env
  Use channels
  Implicit None
  Integer x,y,z,px,py,pz,ch
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: f
  Real(kind=drk), Dimension(1:nx*pex) :: tmp
  Integer :: mpi_err
  Integer npe

#ifdef USE_SIMPLE_IO

! 8/6/2008: M.Shay: Very simple I/O scheme

! Each processor reads in its data from a single file.

  ch=myproc+1
  read (unit(ch)) f(1:nx,1:ny,1:nz)

  call MPI_Barrier(MPI_COMM_WORLD,mpi_err)


#elif USE_OLD_IO

! Old I/O scheme: lines of data are distributed across the output files

  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status

  do pz=0,pez-1; do z=1,nz
    do py=0,pey-1; do y=1,ny
      ch=mod(py*ny+y-1,nchannels)+1
      if (myproc == pe(ch)) then
        read (unit(ch)) tmp
        do px=0,pex-1
          npe = px+pex*py+pex*pey*pz
          if (npe /= myproc) then
             call MPI_Ssend(tmp(px*nx+1),nx,mpi_drk,npe,789,MPI_COMM_WORLD,mpi_err)
          else
             f(1:nx,y,z) = tmp(px*nx+1:px*nx+nx)
          endif
        enddo
      elseif (my_pez==pz .and. my_pey==py) then
        call MPI_Recv(f(1,y,z),nx,mpi_drk,pe(ch),789,MPI_COMM_WORLD,mpi_status,mpi_err)
      endif
    enddo; enddo
  enddo; enddo

#else

! New I/O scheme: slices of data are distributed across the output files

  Integer, Dimension(MPI_STATUS_SIZE,nz) :: mpi_status
  Integer, Dimension(nz) :: mpi_request
  Real(kind=drk), Dimension(0:nx+1,1:ny,0:pex-1) :: tmps

! Every PE sets up a nonblocking receive for all of its slices using MPI_Irecv
! Important:
!   There may be a limit on outstanding MPI_Irecv requests depending
!   on the MPI implementation. If nz is bigger than this limit, the
!   following code may fail !!!!
!   In this case the limit should be increased.
!   (On SGI/Cray systems by setting the environment variable MPI_REQUEST_MAX)

  do z=1,nz
    ch = mod(my_pez*nz+z-1,nchannels)+1
    call MPI_Irecv(f(0,1,z),(nx+2)*ny,mpi_drk,pe(ch),789,MPI_COMM_WORLD,mpi_request(z),mpi_err)
  enddo

! The PEs responsible for input read and scatter the slices

  do ch = 1,nchannels
    if (myproc == pe(ch)) then
      do z=ch,pez*nz,nchannels
        pz = (z-1)/nz
        do py=0,pey-1
          do y=1,ny
            read (unit(ch)) tmp
            do px=0,pex-1
              tmps(1:nx,y,px) = tmp(px*nx+1:px*nx+nx)
            enddo
          enddo
          do px=0,pex-1
            npe = px+pex*py+pex*pey*pz
            call MPI_Send(tmps(0,1,px),(nx+2)*ny,mpi_drk,npe,789,MPI_COMM_WORLD,mpi_err)
          enddo
        enddo
      enddo
    endif
  enddo

! All PEs have to wait until all requests are completed

  call MPI_Waitall(nz,mpi_request,mpi_status,mpi_err)

  call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

#endif


  return
End Subroutine readfield

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifndef twofluid
Subroutine readparticles(rv,np)
  Use pe_env
  Use channels
  Implicit None
  Real(kind=prk), Dimension(6,maxparticles) :: rv
  Real(kind=prk), Dimension(6,bufsize) :: buffer
  Integer :: npp, n_pes_in
  Integer :: np,ppe,pos,ch
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status
  Integer :: mpi_err, dest_pe, num


#ifdef USE_SIMPLE_IO
! 8/6/2008: Very simple I/O. Each processor reads a single file.

  ch=myproc+1
  read(unit(ch)) n_pes_in
  read(unit(ch)) np
  read(unit(ch)) rv(:,1:np)

#else

  np=0
  n_pes_in=0
  do ch=1,nchannels
    if (myproc == pe(ch)) read(unit(ch)) n_pes_in
  enddo

! We use the fact that:  pe(ch)=mod(ch-1,n_pes)

  call MPI_Bcast(n_pes_in,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)

  do ppe=0,n_pes_in-1

    ch=mod(ppe,nchannels)+1

    if (n_pes <= nchannels) then
      dest_pe = pe(ch)
    else
      !
      ! For performance reasons, n_pes should be a multiple of nchannels
      ! in this case, since otherwise a PE which is reading has to send
      ! the data to another PE which is reading, which causes the I/O
      ! to get serialized.
      ! The code should work for any number of n_pes/nchannels, however
      !
      dest_pe = mod(ppe,n_pes)
    endif

    if (myproc == pe(ch)) then
      read(unit(ch)) npp
      if(myproc /= dest_pe) then
        call MPI_Ssend(npp,1,MPI_INTEGER,dest_pe,1234,MPI_COMM_WORLD,mpi_err)
      else
        if(np+npp .gt. maxparticles) then
          write(6,*) '***** input: particle buffer overflow *****'
          call exitallpes()
        endif
      endif
      pos=0
      do while (pos<npp)
        read(unit(ch)) buffer
        num = min(bufsize,npp-pos)
        if(myproc /= dest_pe) then
          call MPI_Ssend(buffer,6*num,mpi_prk,dest_pe,5678,MPI_COMM_WORLD,mpi_err)
        else
          rv(:,np+pos+1:np+pos+num) = buffer(:,:)
        endif
        pos=pos+bufsize
      enddo
      if(myproc == dest_pe) np=np+npp
    elseif(myproc == dest_pe) then
      call MPI_Recv(npp,1,MPI_INTEGER,pe(ch),1234,MPI_COMM_WORLD,mpi_status,mpi_err)
      if(np+npp .gt. maxparticles) then
        write(6,*) '***** input: particle buffer overflow *****'
        call exitallpes()
      endif
      pos=0
      do while (pos<npp)
        num = min(bufsize,npp-pos)
        call MPI_Recv(rv(1,np+pos+1),6*num,mpi_prk,pe(ch),5678,MPI_COMM_WORLD, &
                      mpi_status,mpi_err)
        pos=pos+bufsize
      enddo
      np=np+npp
    endif
  enddo

  call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

#endif

  call redistribute(rv,np)
  return 
End Subroutine readparticles
#endif
