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
End Module channels

!------------------------------------------------------------------------------

Program vdist
  Use channels
  Implicit None
  Character*16 :: filename, argbuf
  Integer :: count,n1,n2,n3,n4,n5,ch,i
  Real(kind=drk), Dimension(-vbin:vbin) :: vx,vy,vz
  Real(kind=drk) :: t,vmaxe,vmaxi,n_avg

  Real(kind=drk),Dimension(0:vbin) :: energy

! read filenumber, vmax and subvolume
  call our_getarg(1,argbuf)
  if (argbuf == ' ') then
    write(6,*) '***** vdist: missing number of input file *****'
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
      write(6,*) '***** vdist: input-file incompatible *****'
! write(6,*) n1,nx,pex,n2,ny,pey,n3,nz,pez,n4,bufsize,n5,nchannels
      stop
    endif
  enddo
! skip field data
  call skipfield()
! read and diagnose ions
  call diagparticles(vx,vy,vz,energy,vmaxi)
  write(6,*)
  write(6,"(A)") '** Ion velocity distribution: v, vx(v), vy(v), vz(v) **'
  do i=-vbin,vbin
    write(6,"(G13.5,E13.5,E13.5,E13.5)") (i*vmaxi)/vbin,vx(i),vy(i),vz(i)
  enddo

! write ion energy
  write(6,*)
  write(6,"(A)") '** Ion energy: Kinetic energy, KE(v) **'
  do i=0,vbin
    write(6,"(G13.5,E13.5)") (i*(vmaxi/sqrt(c_2)))/vbin,energy(i)
  enddo

#ifndef hybrid
! read and diagnose electrons
  call diagparticles(vx,vy,vz,energy,vmaxe)
  write(6,*)
  write(6,"(A)") '** Electron velocity distribution: v, vx(v), vy(v), vz(v) **'
  do i=-vbin,vbin
    write(6,"(G13.5,E13.5,E13.5,E13.5)") (i*vmaxe)/vbin,vx(i),vy(i),vz(i)
  enddo

! write electron energy
  write(6,*)
  write(6,"(A)") '** Electron energy: Kinetic energy, KE(v) **'
  do i=0,vbin
    write(6,"(G13.5,E13.5)") (i*(vmaxe/sqrt(c_2)))/vbin,energy(i)
  enddo

#endif
! close files
  do ch=1,nchannels
    close(unit(ch))
  enddo
End Program vdist

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

Subroutine diagparticles(vx,vy,vz,energy,vmax)
  Use channels
  Implicit None
  Real(kind=drk), Dimension(-vbin:vbin) :: vx,vy,vz
  Real(kind=prk), Dimension(6,maxparticles), Save :: rv
  Real(kind=prk), Dimension(6,bufsize) :: buffer
  Real(kind=drk) :: vmax,tot
  Integer :: np,ppe,pos,n_pes_in,ch,i,bin

  Real(kind=drk), Dimension(0:vbin) :: energy
  Integer :: bin_en
  Real(kind=drk) :: kin_en,vbyc_sq

! delete velocity bins
  do i=-vbin,vbin
    vx(i)=0.; vy(i)=0.; vz(i)=0.
  enddo
  do i = 0,vbin
    energy(i) = 0.
  enddo

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
          write(6,*) '***** vdist: particle buffer overflow *****'
          stop
        endif
    enddo
!       analyse particle buffer
    do i=1,np
      if (rv(1,i) .ge. xmin*lx .and. rv(1,i) .le. xmax*lx .and. &
          rv(2,i) .ge. ymin*ly .and. rv(2,i) .le. ymax*ly .and. &
          rv(3,i) .ge. zmin*lz .and. rv(3,i) .le. zmax*lz) then
        bin=int(abs(rv(4,i))/vmax*vbin+0.5)*sign(real(1.0,prk),rv(4,i))
        if (bin .ge. -vbin .and. bin .le. vbin) vx(bin)=vx(bin)+1.
        bin=int(abs(rv(5,i))/vmax*vbin+0.5)*sign(real(1.0,prk),rv(5,i))
        if (bin .ge. -vbin .and. bin .le. vbin) vy(bin)=vy(bin)+1.
        bin=int(abs(rv(6,i))/vmax*vbin+0.5)*sign(real(1.0,prk),rv(6,i))
        if (bin .ge. -vbin .and. bin .le. vbin) vz(bin)=vz(bin)+1.
! Add bit to calculate energy.  Relativistic KE/mc^2 = (gamma-1)
        vbyc_sq = (rv(4,i)**2+rv(5,i)**2+rv(6,i)**2)/c_2
! Check to see if particle exceeds lightspeed.  Should only happen in 
! non-relativistic simulations.
        if (vbyc_sq .lt. 1.) then
          kin_en = 1./sqrt(1.-vbyc_sq)-1.  
! A hack: Use vmax again as the maximum.  The default (vmax/c^2=1)
! means max(KE) is the rest mass.  To be general there should be another
! flag in the call to vdist, but who wants that?
          bin_en = int(kin_en/(vmax/sqrt(c_2))*vbin+0.5)
          if (bin_en .le. vbin) energy(bin_en)=energy(bin_en)+1.
        endif
      endif
    enddo
  enddo
! normalize distribution functions
! MS: Old version had normalization of any component equal to 
! vbin/vmax = vbin/sqrt(c^2).  Changed to make it 1.
  tot=0.; do i=-vbin,vbin; tot=tot+vx(i); enddo
  if (tot .ne. 0) then
    do i=-vbin,vbin; vx(i)=vx(i)/tot; enddo
  endif
  tot=0.; do i=-vbin,vbin; tot=tot+vy(i); enddo
  if (tot .ne. 0) then
    do i=-vbin,vbin; vy(i)=vy(i)/tot; enddo
  endif
  tot=0.; do i=-vbin,vbin; tot=tot+vz(i); enddo
  if (tot .ne. 0) then
    do i=-vbin,vbin; vz(i)=vz(i)/tot; enddo
  endif  
  tot=0.; do i=0,vbin; tot=tot+energy(i); enddo
  if (tot .ne. 0) then
    do i=0,vbin; energy(i)=energy(i)/tot; enddo
  endif
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
