!******************************************************************************
!                   3D particle code: velocity distribution
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
  Integer :: max_part
End Module channels

!------------------------------------------------------------------------------

Program ener
  Use channels
  Implicit None
  Character*16 :: filename, argbuf
  Integer :: count,n1,n2,n3,n4,n5,ch,i
  Real(kind=drk) :: t,vmaxe,vmaxi,n_avg

  Real(kind=drk) :: e_min
  Real(kind=prk), Allocatable, Dimension(:,:) :: part_arr
  Integer :: number

! read filenumber, vmax and subvolume
  call our_getarg(1,argbuf)
  if (argbuf == ' ') then
    write(6,*) '***** vdist: missing number of input file *****'
    stop
  endif
  read(argbuf,*) count
  e_min=1.; max_part=10; xmin=0.; xmax=1.; ymin=0.; ymax=1.; zmin=0.; zmax=1.
! vmaxe has no meaning in hybrid code, read anyway to keep interface compatible
  call our_getarg(2,argbuf); if (argbuf.ne.' ') read(argbuf,*) e_min
  call our_getarg(3,argbuf); if (argbuf.ne.' ') read(argbuf,*) max_part
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
  write(6,*) 'Minimum energy/mc2:',e_min

! assign input channels and open files
  do ch=1,nchannels
    unit(ch)=30+ch
    write(filename,"(a,'-',i3.3,'.',i3.3)") head, ch, count
    open(unit=unit(ch),file=filename,action='read',form='unformatted')
    read(unit(ch)) t,n_avg,n1,n2,n3,n4,n5
    if (n1 .ne. nx*pex .or. n2 .ne. ny*pey .or. n3 .ne. nz*pez &
                       .or. n4 .ne. bufsize .or. n5 .ne. nchannels) then
      write(6,*) '***** ener: input-file incompatible *****'
! write(6,*) n1,nx,pex,n2,ny,pey,n3,nz,pez,n4,bufsize,n5,nchannels
      stop
    endif
  enddo
! skip field data
  call skipfield()

! Allocate space.  Since max_part is an input I need to make part_arr
! an allocatable array.
  allocate(part_arr(6,max_part))

  number = 1
! read and diagnose ions
  call diagparticles(e_min,part_arr,number)
  write(6,*) number
  do i=1,number-1
    write(6,"(6E13.5)") part_arr(:,i)
  enddo

#ifndef hybrid
! read and diagnose electrons
  call diagparticles(e_min,part_arr,number)
  write(6,*) number
  do i=1,number-1
    write(6,"(6E13.5)") part_arr(:,i)
  enddo
#endif

! close files
  do ch=1,nchannels
    close(unit(ch))
  enddo
End Program ener

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

Subroutine diagparticles(e_min,ener_part,k)
  Use channels
  Implicit None
  Real(kind=prk), Dimension(6,maxparticles), Save :: rv
  Real(kind=prk), Dimension(6,bufsize) :: buffer
  Integer :: np,ppe,pos,n_pes_in,ch,i,bin

  Real(kind=drk) :: kin_en,vbyc_sq
  Integer :: k
  Real(kind=prk), Dimension(6,max_part) :: ener_part
  Real(kind=drk) :: e_min

! Initial array index of energetic particles
  k = 1

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
          write(6,*) '***** ener: particle buffer overflow *****'
          stop
        endif
    enddo
! analyse particle buffer
    if (k .lt. max_part) then 
      do i=1,np
! Add bit to calculate energy.  Relativistic KE/mc^2 = (gamma-1)
        vbyc_sq = (rv(4,i)**2+rv(5,i)**2+rv(6,i)**2)/c_2
! Check to see if particle exceeds lightspeed.  Should only happen in 
! non-relativistic simulations.
        if (vbyc_sq .lt. 1.) then
          kin_en = 1./sqrt(1.-vbyc_sq)-1.  
          if (rv(1,i) .ge. xmin*lx .and. rv(1,i) .le. xmax*lx .and. &
            rv(2,i) .ge. ymin*ly .and. rv(2,i) .le. ymax*ly .and. &
            rv(3,i) .ge. zmin*lz .and. rv(3,i) .le. zmax*lz .and. &
            kin_en .ge. e_min) then
              ener_part(:,k) = rv(:,i)
              k = k+1
	      if (k .gt. max_part) exit
          endif
        endif
      enddo
    endif
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
