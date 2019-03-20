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
!					Output distribution functions
!------------------------------------------------------------------------------

Subroutine output_dist(fdist,nvx,fh,unit2,offset)

  Use pe_env
  Use timer
  Implicit None

  Integer :: unit2,fh,x,y,z,nvx,vx
  Real(kind=drk), Dimension(nx,nvx) :: fdist
  Integer(kind=MPI_OFFSET_KIND) :: offset

  Character, Dimension(2*nx,nvx) :: char_fdist
  Real(kind=drk) :: lmin, lmax
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status
  Integer :: mpi_err,digit0,digit1,intdbl

!  call start_timer()
! calculate global minimum and maximum
  lmin = 1e20; lmax = -1e20
  do vx=1,nvx ; do x=1,nx
    lmin=min(lmin,fdist(x,vx)); lmax=max(lmax,fdist(x,vx))
  enddo; enddo
  lmax=totalmax(lmax); lmin=totalmin(lmin)

  if (abs(lmin-lmax) .lt. eps) then
    char_fdist = char(0)
  else
    do vx=1,nvx ; do x=1,nx
       intdbl = nint((fdist(x,vx)-lmin)/(lmax-lmin)*65535)
       digit0 = iand(intdbl,255)  ! First 8 bytes
       digit1 = ishft(iand(intdbl,65280),-8) ! Next 8 bytes
       char_fdist(2*x-1,vx) = char(digit1) ! like a readl 16 bit output
       char_fdist(2*x,vx) = char(digit0) 

    enddo; enddo
  endif
  if (myproc == 0) write(unit2,"(E14.6,E14.6)") lmin,lmax

  call MPI_FILE_WRITE_AT_ALL(fh,offset,char_fdist,2*nx*nvx,MPI_CHARACTER,&
                             mpi_status,mpi_err)

! Not necessary according to Johnathan Carter at NERSC
!  call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

!  call stop_timer(t_movie,n_movie)
  return
End Subroutine output_dist
