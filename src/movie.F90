!******************************************************************************
!                      3D particle code: movie output
!                           Andreas Zeiler, 1999
!
!                        LATEST CHANGE: May 3, 2001
!
!******************************************************************************

#include "param"
#ifndef double_byte
#define MIN 0
#define STEPS 255
#else
#define MIN (-32768)
#define STEPS 65535
#endif

!------------------------------------------------------------------------------
!        				full resolution movie: x,y planes 
!------------------------------------------------------------------------------

Subroutine moviexyfull(f,z,unit1,frmtd,line)
! write xy plane with given z
  Use pe_env
  Use timer
  Implicit None

  Integer x,y,z,unit1,px,py,pz,pe,zz,yy,i,j
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: f
#ifdef double_precision
  Real(kind=drk), Dimension(nx,ny) :: plane2d
  Real(kind=drk), Dimension(nx) :: tmp
#else
  Real(kind=prk), Dimension(nx,ny) :: plane2d
  Real(kind=prk), Dimension(nx) :: tmp
#endif
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status
  Integer :: mpi_err
  Logical :: frmtd,line

  call start_timer()
  pz=(z-1)/nz; zz=mod(z-1,nz)+1
      do y=1,ny; do x=1,nx
        plane2d(x,y)=f(x,y,zz)
      enddo; enddo
! output data

 if (line) then
    y = ny*(pey/2 + 1 )
    py=(y-1)/ny; yy=mod(y-1,ny)+1
    do px=0,pex-1
      pe=px+py*pex+pz*pex*pey
      if (myproc == 0) then
        if (pe /= 0) then
          !call MPI_Recv(tmp,nx,MPI_INTEGER,pe,11,MPI_COMM_WORLD,mpi_status,mpi_err)
#ifdef double_precision
          call MPI_Recv(tmp,nx,MPI_drk,pe,11,MPI_COMM_WORLD,mpi_status,mpi_err)
#else
          call MPI_Recv(tmp,nx,MPI_prk,pe,11,MPI_COMM_WORLD,mpi_status,mpi_err)
#endif
        else
          tmp(:) = plane2d(:,yy)
        endif

        if (frmtd) then
         write(unit1,FMT="(E14.6)") tmp
        else
         call writelinefull(unit1,tmp,nx)
        endif

      elseif (myproc == pe) then
        !call MPI_Ssend(plane2d(1,yy),nx,MPI_INTEGER,0,11,MPI_COMM_WORLD,mpi_err)
#ifdef double_precision
        call MPI_Ssend(plane2d(1,yy),nx,MPI_drk,0,11,MPI_COMM_WORLD,mpi_err)
#else
        call MPI_Ssend(plane2d(1,yy),nx,MPI_prk,0,11,MPI_COMM_WORLD,mpi_err)
#endif
      endif
    enddo

 else

  do y=1,ny*pey
    py=(y-1)/ny; yy=mod(y-1,ny)+1
    do px=0,pex-1
      pe=px+py*pex+pz*pex*pey
      if (myproc == 0) then
        if (pe /= 0) then
          !call MPI_Recv(tmp,nx,MPI_INTEGER,pe,11,MPI_COMM_WORLD,mpi_status,mpi_err)
#ifdef double_precision
          call MPI_Recv(tmp,nx,MPI_drk,pe,11,MPI_COMM_WORLD,mpi_status,mpi_err)
#else
          call MPI_Recv(tmp,nx,MPI_prk,pe,11,MPI_COMM_WORLD,mpi_status,mpi_err)
#endif
        else
          tmp(:) = plane2d(:,yy)
        endif

        if (frmtd) then
         write(unit1,FMT="(E14.6)") tmp
        else
         call writelinefull(unit1,tmp,nx)
        endif

      elseif (myproc == pe) then
        !call MPI_Ssend(plane2d(1,yy),nx,MPI_INTEGER,0,11,MPI_COMM_WORLD,mpi_err)
#ifdef double_precision
        call MPI_Ssend(plane2d(1,yy),nx,MPI_drk,0,11,MPI_COMM_WORLD,mpi_err)
#else
        call MPI_Ssend(plane2d(1,yy),nx,MPI_prk,0,11,MPI_COMM_WORLD,mpi_err)
#endif
      endif
    enddo
  enddo

  endif ! End the IF statement for line output.

  call MPI_Barrier(MPI_COMM_WORLD,mpi_err) ! Not absolutly necessary
  call stop_timer(t_movie,n_movie)
  return
End Subroutine moviexyfull

!------------------------------------------------------------------------------
!							  movie: x,y planes 
!------------------------------------------------------------------------------

Subroutine moviexy(f,z,unit1,unit2)
! write xy plane with given z
  Use pe_env
  Use timer
  Implicit None

  Integer x,y,z,unit1,unit2,px,py,pz,pe,zz,yy
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: f
  Integer, Dimension(nx,ny) :: plane2d
  Integer, Dimension(nx) :: tmp
  Real(kind=drk) :: lmin, lmax
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status
  Integer :: mpi_err

  call start_timer()
  pz=(z-1)/nz; zz=mod(z-1,nz)+1
! calculate global minimum and maximum
  lmin = 1e20; lmax = -1e20
  if (pz == my_pez) then
    do y=1,ny; do x=1,nx
      lmin=min(lmin,f(x,y,zz)); lmax=max(lmax,f(x,y,zz))
    enddo; enddo
  endif
  lmax=totalmax(lmax); lmin=totalmin(lmin)
! byte scale data
  if (pz == my_pez) then
    if (abs(lmin-lmax) .lt. eps) then
      plane2d=0.0
    else
      do y=1,ny; do x=1,nx
        plane2d(x,y)=nint((f(x,y,zz)-lmin)/(lmax-lmin)*STEPS)+MIN
      enddo; enddo
    endif
  endif
! output data
  call MPI_Barrier(MPI_COMM_WORLD,mpi_err) ! Not absolutly necessary
  if (myproc == 0) then
    write(unit2,"(E14.6,E14.6)") lmin,lmax

! An XLF routine that flushes the I/O buffer. This way on a code crash
! the log file still ha1s something useful in it.
#   ifdef AIX
      call flush_(unit2)
#   endif
#   ifdef Linux
      call flush(unit2)
#   endif

  endif
  do y=1,ny*pey
    py=(y-1)/ny; yy=mod(y-1,ny)+1
    do px=0,pex-1
      pe=px+py*pex+pz*pex*pey
      if (myproc == 0) then
        if (pe /= 0) then
          call MPI_Recv(tmp,nx,MPI_INTEGER,pe,11,MPI_COMM_WORLD,mpi_status,mpi_err)
        else
          tmp(:) = plane2d(:,yy)
        endif
        call writeline(unit1,tmp,nx)
      elseif (myproc == pe) then
        call MPI_Ssend(plane2d(1,yy),nx,MPI_INTEGER,0,11,MPI_COMM_WORLD,mpi_err)
      endif
    enddo
  enddo
  call MPI_Barrier(MPI_COMM_WORLD,mpi_err) ! Not absolutly necessary
  call stop_timer(t_movie,n_movie)
  return
End Subroutine moviexy

!------------------------------------------------------------------------------
!							  movie: x,z planes 
!------------------------------------------------------------------------------

Subroutine moviexz(f,y,unit1,unit2)
! write xz plane with given y
  Use pe_env
  Use timer
  Implicit None

  Integer x,y,z,unit1,unit2,px,py,pz,pe,yy,zz
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: f
  Integer, Dimension(nx,nz) :: plane2d
  Integer, Dimension(nx) :: tmp
  Real(kind=drk) :: lmin, lmax
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status
  Integer :: mpi_err

  call start_timer()
  py=(y-1)/ny; yy=mod(y-1,ny)+1
! calculate global minimum and maximum
  lmin = 1e20; lmax = -1e20
  if (py == my_pey) then
    do z=1,nz; do x=1,nx
      lmin=min(lmin,f(x,yy,z)); lmax=max(lmax,f(x,yy,z))
    enddo; enddo
  endif
  lmax=totalmax(lmax); lmin=totalmin(lmin)
! byte scale data
  if (py == my_pey) then
    if (abs(lmin-lmax) .lt. eps) then
      plane2d=0.0
    else
      do z=1,nz; do x=1,nx
        plane2d(x,z)=nint((f(x,yy,z)-lmin)/(lmax-lmin)*STEPS)+MIN
      enddo; enddo
    endif
  endif
! output data
  call MPI_Barrier(MPI_COMM_WORLD,mpi_err) ! Not absolutly necessary
  if (myproc == 0) write(unit2,"(E14.6,E14.6)") lmin,lmax
  do z=1,nz*pez
    pz=(z-1)/nz; zz=mod(z-1,nz)+1
    do px=0,pex-1
      pe=px+py*pex+pz*pex*pey
      if (myproc == 0) then
        if (pe /= 0) then
          call MPI_Recv(tmp,nx,MPI_INTEGER,pe,22,MPI_COMM_WORLD,mpi_status,mpi_err)
        else
          tmp(:) = plane2d(:,zz)
        endif
        call writeline(unit1,tmp,nx)
      elseif (myproc == pe) then
        call MPI_Ssend(plane2d(1,zz),nx,MPI_INTEGER,0,22,MPI_COMM_WORLD,mpi_err)
      endif
    enddo
  enddo
  call MPI_Barrier(MPI_COMM_WORLD,mpi_err) ! Not absolutly necessary
  call stop_timer(t_movie,n_movie)
  return
End Subroutine moviexz

!------------------------------------------------------------------------------
!							  movie: y,z planes 
!------------------------------------------------------------------------------

Subroutine movieyz(f,x,unit1,unit2)
! write yz plane with given x
  Use pe_env
  Use timer
  Implicit None

  Integer x,y,z,unit1,unit2,px,py,pz,pe,xx,zz
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: f
  Integer, Dimension(ny,nz) :: plane2d
  Integer, Dimension(ny) :: tmp
  Real(kind=drk) :: lmin, lmax
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status
  Integer :: mpi_err

  call start_timer()
  px=(x-1)/nx; xx=mod(x-1,nx)+1
! calculate global minimum and maximum
  lmin = 1e20; lmax = -1e20
  if (px == my_pex) then
    do z=1,nz; do y=1,ny
      lmin=min(lmin,f(xx,y,z)); lmax=max(lmax,f(xx,y,z))
    enddo; enddo
  endif
  lmax=totalmax(lmax); lmin=totalmin(lmin)
! byte scale data
  if (px == my_pex) then
    if (abs(lmin-lmax) .lt. eps) then
      plane2d=0.0
    else
      do z=1,nz; do y=1,ny
        plane2d(y,z)=nint((f(xx,y,z)-lmin)/(lmax-lmin)*STEPS)+MIN
      enddo; enddo
    endif
  endif
! output data
  call MPI_Barrier(MPI_COMM_WORLD,mpi_err) ! Not absolutly necessary
  if (myproc == 0) write(unit2,"(E14.6,E14.6)") lmin,lmax
  do z=1,nz*pez
    pz=(z-1)/nz; zz=mod(z-1,nz)+1
    do py=0,pey-1
      pe=px+py*pex+pz*pex*pey
      if (myproc == 0) then
        if (pe /= 0) then
          call MPI_Recv(tmp,ny,MPI_INTEGER,pe,33,MPI_COMM_WORLD,mpi_status,mpi_err)
        else
          tmp(:) = plane2d(:,zz)
        endif
        call writeline(unit1,tmp,ny)
      elseif (myproc == pe) then
        call MPI_Ssend(plane2d(1,zz),ny,MPI_INTEGER,0,33,MPI_COMM_WORLD,mpi_err)
      endif
    enddo
  enddo
  call MPI_Barrier(MPI_COMM_WORLD,mpi_err) ! Not absolutly necessary
  call stop_timer(t_movie,n_movie)
  return
End Subroutine movieyz



!------------------------------------------------------------------------------
!                                                         volume: x,y,z volume
!------------------------------------------------------------------------------

Subroutine volume(f,unit1,unit2)
! write xyz volume
  Use pe_env
  Use timer
  Implicit None

  Integer x,y,z,unit1,unit2,px,py,pz,pe,zz,yy
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: f
  Integer, Dimension(nx,ny) :: plane2d
  Integer, Dimension(nx) :: tmp
  Real(kind=drk) :: lmin, lmax
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status
  Integer :: mpi_err

  call start_timer()

! calculate global minimum and maximum
  lmin = 1e20; lmax = -1e20
  do z=1,nz; do y=1,ny; do x=1,nx
    lmin=min(lmin,f(x,y,z)); lmax=max(lmax,f(x,y,z))
  enddo; enddo; enddo
  lmax=totalmax(lmax); lmin=totalmin(lmin)
  if (myproc == 0) write(unit2,"(E14.6,E14.6)") lmin,lmax
! loop over z planes
  do z=1,nz*pez
    pz=(z-1)/nz; zz=mod(z-1,nz)+1
    if (pz == my_pez) then
      if (abs(lmin-lmax) .lt. eps) then
        plane2d=0.0
      else
        do y=1,ny; do x=1,nx
          plane2d(x,y)=nint((f(x,y,zz)-lmin)/(lmax-lmin)*STEPS)+MIN
        enddo; enddo
      endif
    endif
    call MPI_Barrier(MPI_COMM_WORLD,mpi_err) ! Not absolutly necessary
    do y=1,ny*pey
      py=(y-1)/ny; yy=mod(y-1,ny)+1
      do px=0,pex-1
        pe=px+py*pex+pz*pex*pey
        if (myproc == 0) then
          if (pe /= 0) then
            call MPI_Recv(tmp,nx,MPI_INTEGER,pe,44,MPI_COMM_WORLD,mpi_status,mpi_err)
          else
            tmp(:) = plane2d(:,yy)
          endif
          call writeline(unit1,tmp,nx)
        elseif (myproc == pe) then
          call MPI_Ssend(plane2d(1,yy),nx,MPI_INTEGER,0,44,MPI_COMM_WORLD,mpi_err)
        endif
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD,mpi_err) ! Not absolutly necessary
  enddo
  call stop_timer(t_movie,n_movie)
  return
End Subroutine volume


!------------------------------------------------------------------------------
!					volume_z: x,y,z volume of one z level
!------------------------------------------------------------------------------

Subroutine volume_z(f,fh,unit2,offset)
! write xyz volume of one zlevel of processors
  Use pe_env
  Use timer
  Implicit None

  Integer :: unit2,fh,x,y,z
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: f
  Integer(kind=MPI_OFFSET_KIND) :: offset

  Character, Dimension(nx,ny,nz) :: char_f
  Real(kind=drk) :: lmin, lmax
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status
  Integer :: mpi_err

  call start_timer()
! calculate global minimum and maximum
  lmin = 1e20; lmax = -1e20
  do z=1,nz; do y=1,ny; do x=1,nx
    lmin=min(lmin,f(x,y,z)); lmax=max(lmax,f(x,y,z))
  enddo; enddo; enddo
  lmax=totalmax(lmax); lmin=totalmin(lmin)

  if (abs(lmin-lmax) .lt. eps) then
    char_f = char(0)
  else
    do z = 1,nz; do y = 1,ny; do x = 1,nx
      char_f(x,y,z)= char(nint((f(x,y,z)-lmin)/(lmax-lmin)*STEPS)+MIN)
    enddo; enddo; enddo
  endif
  if (myproc == 0) write(unit2,"(E14.6,E14.6)") lmin,lmax

  call MPI_FILE_WRITE_AT_ALL(fh,offset,char_f,nx*ny*nz,MPI_CHARACTER,&
                             mpi_status,mpi_err)

! Not necessary according to Johnathan Carter at NERSC
!  call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

  call stop_timer(t_movie,n_movie)
  return
End Subroutine volume_z

!------------------------------------------------------------------------------
!                   Simple movie output on each processor
!------------------------------------------------------------------------------
Subroutine moviesimple(f,fh,unit2,recnum)
! Each processor writes its own file
  Use pe_env
  Use timer
  Implicit None

  Integer :: unit2,fh,x,y,z
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: f
  Integer :: recnum

  Character, Dimension(nx,ny,nz) :: char_f
  Real(kind=drk) :: lmin, lmax
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status
  Integer :: mpi_err

  call start_timer()
! calculate global minimum and maximum
  lmin = 1e20; lmax = -1e20
  do z=1,nz; do y=1,ny; do x=1,nx
    lmin=min(lmin,f(x,y,z)); lmax=max(lmax,f(x,y,z))
  enddo; enddo; enddo
  lmax=totalmax(lmax); lmin=totalmin(lmin)

  if (abs(lmin-lmax) .lt. eps) then
    char_f = char(0)
  else
    do z = 1,nz; do y = 1,ny; do x = 1,nx
      char_f(x,y,z)= char(nint((f(x,y,z)-lmin)/(lmax-lmin)*STEPS)+MIN)
    enddo; enddo; enddo
  endif
  if (myproc == 0) write(unit2,"(E14.6,E14.6)") lmin,lmax

  write(unit=fh,rec=recnum) char_f
  recnum=recnum+1

  call stop_timer(t_movie,n_movie)
  return
End Subroutine moviesimple

!------------------------------------------------------------------------------
!                   Simple movie output on each processor
!------------------------------------------------------------------------------
Subroutine moviesimplefull(f,fh,recnum)
! Each processor writes its own file
  Use pe_env
  Use timer
  Implicit None

  Integer :: fh
  Real(kind=drk), Dimension(0:nx+1,0:ny+1,0:nz+1) :: f
  Integer :: recnum

  Real(kind=drk), Dimension(nx,ny,nz) :: f1
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status
  Integer :: mpi_err

  call start_timer()
  f1=f(1:nx,1:ny,1:nz)
  write(unit=fh,rec=recnum) f1
  recnum=recnum+1

  call stop_timer(t_movie,n_movie)
  return
End Subroutine moviesimplefull

!------------------------------------------------------------------------------
!                          support routine: writeline
!------------------------------------------------------------------------------

Module movie_io
  Integer, Dimension(99) :: rec_ctr
End Module movie_io

Subroutine writeline(unit,tmp,nn)

  Use movie_io

  Implicit None

  Integer :: x,unit,outh,outl,nn
  Integer, Dimension(nn) :: tmp
#ifdef double_byte
  Character, Dimension(2*nn) :: ctmp
#else
  Character, Dimension(nn) :: ctmp
#endif

  do x=1,nn
#ifndef double_byte
    ctmp(x) = char(tmp(x))
#else
    outh=abs(tmp(x)/256) 
    outl=mod(tmp(x),256)
    if (outl<0) outh=255-outh
    if (outl==0 .and. tmp(x)<0) outh=256-outh
    ctmp(2*x-1) = char(outh)
    ctmp(2*x  ) = char(outl)
#endif
  enddo

  rec_ctr(unit) = rec_ctr(unit)+1
  write(unit,rec=rec_ctr(unit)) ctmp

  return
End Subroutine writeline

Subroutine reset_rec_ctr
  Use movie_io
  Implicit None

   rec_ctr = 0

End Subroutine reset_rec_ctr



Subroutine writelinefull(unit,tmp,nn)

  Use movie_io

  Implicit None

  Integer :: unit,nn
#ifdef double_precision
  Real(kind=drk), Dimension(nn) :: tmp
#else
  Real(kind=prk), Dimension(nn) :: tmp
#endif
  rec_ctr(unit) = rec_ctr(unit)+1
  write(unit,rec=rec_ctr(unit)) tmp

  return
End Subroutine writelinefull





