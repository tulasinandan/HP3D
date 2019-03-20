! 8/28/2008: M. Shay: Program to combine different processor files
!            when each processor writes its own movie file. 

#include "param"


Program moviecombine3

Implicit None

Integer :: ival, iproc, offset, recin, recout, ntimes, istart, istop, iskip
Integer :: my_pex, my_pey, my_pez, x,y,z, itime, stepnum
Integer, parameter :: nprocs=pex*pey*pez
Character(len=16) :: filename, buf
Integer, Dimension(nprocs) :: mypexA, mypeyA, mypezA, iminA, imaxA &
               ,jminA, jmaxA, kminA, kmaxA
#define moviecombine_variables_full
#include movie_header_full
#undef moviecombine_variables_full
Real(kind=8), Dimension(nx,ny,nz,nvalues) :: chardata
Real(kind=8), Dimension(nx*pex,ny*pey,nz*pez,nvalues) :: charfull

#define movie_dataout_array_full
#include movie_header_full
#undef movie_dataout_array_full

! Get the number of movie files form the command line
! Also get start stop and skip for movie times
!call getarg(1,buf)  !call getarg(1,buf)
!read (buf,*) stepnum!read (buf,*) stepnum
 call getarg(1,buf)  !call getarg(2,buf)
 read (buf,*) ntimes !read (buf,*) ntimes
 call getarg(2,buf)  !call getarg(3,buf)
 read (buf,*) istart !read (buf,*) istart
 call getarg(3,buf)  !call getarg(4,buf)
 read (buf,*) istop  !read (buf,*) istop
 call getarg(4,buf)  !call getarg(5,buf)
 read (buf,*) iskip  !read (buf,*) iskip

print *,'ntimes = ',ntimes
print *,'istart, istop, iskip = ',istart,istop,iskip

!Open destination movie files

do ival=1,nvalues
  open(unit=10+ival,file=valoutputfull(ival),status='unknown',form='unformatted' &
      ,access='direct',recl=8*nx*pex*ny*pey*nz*pez)
end do

! Open raw movie files to read in
do iproc=0,nprocs-1
! write(filename,"(a,'.',i4.4,'.',i1.1)") 'moviefull',iproc,stepnum
  write(filename,"(a,'.',i6.6)") 'moviefull',iproc
  open (unit=100+iproc,file=filename,status='old',action='read',form='unformatted' &
       ,access='direct',recl=8*nx*ny*nz*nvalues)
end do

! Initialize data information and structure for each processor

do iproc=0,nprocs-1
    mypexA(iproc) = Mod(iproc,pex)
    mypeyA(iproc) = Mod(iproc,pex*pey)/pex
    mypezA(iproc) = iproc/(pex*pey)
    iminA(iproc) = mypexA(iproc)*nx+1
    imaxA(iproc) = (mypexA(iproc)+1)*nx
    jminA(iproc) = mypeyA(iproc)*ny+1
    jmaxA(iproc) = (mypeyA(iproc)+1)*ny
    kminA(iproc) = mypezA(iproc)*nz+1
    kmaxA(iproc) = (mypezA(iproc)+1)*nz
end do


!Read in the data and output it

recout=1
do itime=istart,istop,iskip
  do iproc=0,nprocs-1
    !Read in all data for one time
    read(unit=100+iproc,rec=itime) chardata
    !Redistribute data over full grid
    charfull(iminA(iproc):imaxA(iproc),jminA(iproc):jmaxA(iproc),kminA(iproc):kmaxA(iproc),:)&
         = chardata
    if (mod(iproc,100) .eq. 0)  &
       print *,'istep = ',stepnum,'  itime = ',itime,'  iproc = ',iproc
  end do !iproc

  print *,'Writing itime = ',itime
  do ival=1,nvalues

    write(unit=10+ival,rec=recout) charfull(:,:,:,ival)
  end do !ival

  recout=recout+1

end do !itime

do ival=1,nvalues
  close(10+ival)
end do

do iproc=0,nprocs-1
  close(100+iproc)
end do

End Program moviecombine3

