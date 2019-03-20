!******************************************************************************
!               3D particle code: 2-D velocity distribution
!                           Andreas Zeiler, 1999
!
!                        LATEST CHANGE: May 3, 2001
!
!******************************************************************************

#include "param"
#define vbin 100

Module channels
  Implicit  None
  Integer, Dimension(nchannels) :: unit
  Real(kind=drk) :: xmin,xmax,ymin,ymax,zmin,zmax
End Module channels

!------------------------------------------------------------------------------

Program vdist2d
  Use channels
  Implicit None
  Character*16 :: filename, argbuf
  Integer :: count,n1,n2,n3,n4,n5,ch,i,j
  Real(kind=drk), Dimension(-vbin:vbin,-vbin:vbin) :: vxy,vyz,vxz
  Real(kind=drk) :: t,vmaxe,vmaxi,n_avg

! read filenumber, vmax and subvolume
  call our_getarg(1,argbuf)
  if (argbuf == ' ') then
    write(6,*) '***** vdist2d: missing number of input file *****'
    stop
  endif
  read(argbuf,*) count
  vmaxi=1.; vmaxe=1.; xmin=0.; xmax=1.; ymin=0.; ymax=1.; zmin=0.; zmax=1.
! vmaxe has no meaning in hybrid code, read anyway to keep interface compatible
  call our_getarg(2,argbuf); if (argbuf.ne.' ') read(argbuf,*) vmaxi
  call our_getarg(3,argbuf); if (argbuf.ne.' ') read(argbuf,*) vmaxe
  call our_getarg(4,argbuf); if (argbuf.ne.' ') read(argbuf,*) xmin
  call our_getarg(5,argbuf); if (argbuf.ne.' ') read(argbuf,*) xmax
  call our_getarg(6,argbuf); if (argbuf.ne.' ') read(argbuf,*) ymin
  call our_getarg(7,argbuf); if (argbuf.ne.' ') read(argbuf,*) ymax
  call our_getarg(8,argbuf); if (argbuf.ne.' ') read(argbuf,*) zmin
  call our_getarg(9,argbuf); if (argbuf.ne.' ') read(argbuf,*) zmax
  write(filename,"(a,'-xxx.',i3.3)") head, count
  write(6,*) 'input from ', filename
  write(6,*) 'diagnosed volume:'
  write(6,*) '    ',xmin,'<= x/lx <=',xmax
  write(6,*) '    ',ymin,'<= y/ly <=',ymax
  write(6,*) '    ',zmin,'<= z/lz <=',zmax
#ifndef hybrid
  write(6,*) 'vmax_i/c:',vmaxi
  write(6,*) 'vmax_e/c:',vmaxe
  vmaxi=vmaxi*sqrt(c_2); vmaxe=vmaxe*sqrt(c_2)
#else
  write(6,*) 'vmax_i:',vmaxi
#endif
! assign input channels and open files
  do ch=1,nchannels
    unit(ch)=30+ch
    write(filename,"(a,'-',i3.3,'.',i3.3)") head, ch, count
    open(unit=unit(ch),file=filename,action='read',form='unformatted')
    read(unit(ch)) t,n_avg,n1,n2,n3,n4,n5
    if (n1 .ne. nx*pex .or. n2 .ne. ny*pey .or. n3 .ne. nz*pez &
                       .or. n4 .ne. bufsize .or. n5 .ne. nchannels) then
      write(6,*) '***** vdist2d: input-file incompatible *****'
      stop
    endif
  enddo
! skip field data
  call skipfield()
! read and diagnose ions
  call diagparticles(vxy,vyz,vxz,vmaxi)
! write distribution
  write(filename,"('vdixy.',i3.3)") count
  open(unit=20,file=filename,action='write')
  write(filename,"('vdiyz.',i3.3)") count
  open(unit=21,file=filename,action='write')
  write(filename,"('vdixz.',i3.3)") count
  open(unit=22,file=filename,action='write')
  do j=-vbin,vbin; do i=-vbin,vbin
    write(20,*) vxy(i,j); write(21,*) vyz(i,j); write(22,*) vxz(i,j)
  enddo; enddo
  close(20); close(21); close(22)
#ifndef hybrid
! read and diagnose electrons
  call diagparticles(vxy,vyz,vxz,vmaxe)
! write distribution
  write(filename,"('vdexy.',i3.3)") count
  open(unit=20,file=filename,action='write')
  write(filename,"('vdeyz.',i3.3)") count
  open(unit=21,file=filename,action='write')
  write(filename,"('vdexz.',i3.3)") count
  open(unit=22,file=filename,action='write')
  do j=-vbin,vbin; do i=-vbin,vbin
    write(20,*) vxy(i,j); write(21,*) vyz(i,j); write(22,*) vxz(i,j)
  enddo; enddo
  close(20); close(21); close(22)
#endif
! close files
  do ch=1,nchannels
    close(unit(ch))
  enddo
End Program vdist2d

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
      read (unit(ch)) tmp
    enddo; enddo; enddo; enddo
  enddo
  return
End Subroutine skipfield

!------------------------------------------------------------------------------

Subroutine diagparticles(vxy,vyz,vxz,vmax)
  Use channels
  Implicit None
  Real(kind=drk), Dimension(-vbin:vbin,-vbin:vbin) :: vxy,vyz,vxz
  Real(kind=prk), Dimension(6,maxparticles),Save :: rv
  Real(kind=prk), Dimension(6,bufsize) :: buffer
  Real(kind=drk) :: vmax,tot
  Integer :: np,ppe,pos,n_pes_in,ch,i,j,binx,biny,binz

! delete velocity bins
  do i=-vbin,vbin; do j=-vbin,vbin
    vxy(i,j)=0.; vyz(i,j)=0.; vxz(i,j)=0.
  enddo; enddo
! read number of original PEs
  do ch=1,nchannels
    read(unit(ch)) n_pes_in
  enddo
! for each original PE read complete particle buffer and analyse
  do ppe=0,n_pes_in-1
!       read particle buffer of one PE of original run
    ch=mod(ppe,nchannels)+1
    read(unit(ch)) np
    pos=0
    do while (pos<np)
      read(unit(ch)) buffer
      rv(:,pos+1:pos+bufsize) = buffer(:,:)
      pos=pos+bufsize
        if (pos+bufsize .gt. maxparticles) then
          write(6,*) '***** vdist2d: particle buffer overflow *****'
          stop
        endif
    enddo
!       analyse particle buffer
    do i=1,np
      if (rv(1,i) .ge. xmin*lx .and. rv(1,i) .le. xmax*lx .and. &
          rv(2,i) .ge. ymin*ly .and. rv(2,i) .le. ymax*ly .and. &
          rv(3,i) .ge. zmin*lz .and. rv(3,i) .le. zmax*lz) then
        binx=int(abs(rv(4,i))/vmax*vbin+0.5)*sign(real(1.0,prk),rv(4,i))
        biny=int(abs(rv(5,i))/vmax*vbin+0.5)*sign(real(1.0,prk),rv(5,i))
        binz=int(abs(rv(6,i))/vmax*vbin+0.5)*sign(real(1.0,prk),rv(6,i))
        if (abs(binx) .le. vbin .and. abs(biny) .le. vbin) &
           vxy(binx,biny)=vxy(binx,biny)+1.
        if (abs(biny) .le. vbin .and. abs(binz) .le. vbin) &
           vyz(biny,binz)=vyz(biny,binz)+1.
        if (abs(binx) .le. vbin .and. abs(binz) .le. vbin) &
           vxz(binx,binz)=vxz(binx,binz)+1.
      endif
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
#ifdef Linux
  call getarg(num,buf)
#endif

End
