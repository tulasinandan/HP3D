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

Program vdist_par
  Use pe_env
  Use channels
  Implicit None
  Character*16 :: filename, argbuf
  Integer :: count,n1,n2,n3,n4,n5,ch,i
  Real(kind=drk), Dimension(-vbin:vbin) :: vix,viy,viz,vex,vey,vez
  Real(kind=drk) :: t,vmaxe,vmaxi,n_avg

  Real(kind=drk),Dimension(0:vbin) :: ienergy,eenergy,sumien,sumeen
  Integer :: mpi_err,loopcount,tempch
  Real(kind=drk), Dimension(-vbin:vbin) :: sumvix,sumviy,sumviz,sumvex,sumvey,sumvez
  Real(kind=drk) :: tot

  call init_pe_env()

! read filenumber, vmax and subvolume
  call our_getarg(1,argbuf)
  if (argbuf == ' ') then
    write(6,*) '***** vdist_par: missing number of input file *****'
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

  if (myproc == 0) then
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
  endif

  call MPI_Bcast(vmaxi,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
#ifndef hybrid
  call MPI_Bcast(vmaxe,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
#endif
  call MPI_Bcast(xmin,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
  call MPI_Bcast(xmax,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
  call MPI_Bcast(ymin,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
  call MPI_Bcast(ymax,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
  call MPI_Bcast(zmin,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
  call MPI_Bcast(zmax,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)

  sumvix = 0.; sumviy = 0.; sumviz = 0.; sumvex = 0.; sumvey = 0. 
  sumvez = 0.; sumeen = 0.; sumien = 0.
  do loopcount = 1,ceiling(nchannels*1./nprocs)
    ch = (myproc+1) + (loopcount-1)*nprocs
    if (ch .gt. nchannels) exit 
    unit(ch)=30+ch
    write(filename,"(a,'-',i3.3,'.',i3.3)") head, ch, count
    open(unit=unit(ch),file=filename,action='read',form='unformatted')
    read(unit(ch)) t,n_avg,n1,n2,n3,n4,n5
    if (n1 .ne. nx*pex .or. n2 .ne. ny*pey .or. n3 .ne. nz*pez &
                     .or. n4 .ne. bufsize .or. n5 .ne. nchannels) then
      write(6,*) '***** vdist_par: input-file incompatible *****'
 write(6,*) n1,nx,pex,n2,ny,pey,n3,nz,pez,n4,bufsize,n5,nchannels
      call MPI_Finalize(mpi_err)
      stop
    endif

! skip field data
    if (ch == 1) then
      call skipfield()
    endif

! read and diagnose ions
    call diagparticles(vix,viy,viz,ienergy,vmaxi,ch)
    sumvix = sumvix+vix; sumviy = sumviy+viy; sumviz = sumviz+viz
    sumien = sumien +ienergy

#ifndef hybrid
! read and diagnose electrons
    call diagparticles(vex,vey,vez,eenergy,vmaxe,ch)
    sumvex = sumvex+vex; sumvey = sumvey+vey; sumvez = sumvez+vez
    sumeen = sumeen+eenergy
#endif

    close(unit(ch))
  enddo

! Combine data. Note that I get tricky here by reusing vix, viy, vex, etc.
! Here they become the _result_ of summing sumvix over all processors.
  call MPI_Allreduce(sumvix,vix,2*vbin+1,mpi_drk,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(sumvex,vex,2*vbin+1,mpi_drk,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(sumviy,viy,2*vbin+1,mpi_drk,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(sumvey,vey,2*vbin+1,mpi_drk,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(sumviz,viz,2*vbin+1,mpi_drk,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(sumvez,vez,2*vbin+1,mpi_drk,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(sumien,ienergy,vbin+1,mpi_drk,MPI_SUM,MPI_COMM_WORLD,mpi_err)
  call MPI_Allreduce(sumeen,eenergy,vbin+1,mpi_drk,MPI_SUM,MPI_COMM_WORLD,mpi_err)

! Normalize distribution functions
  if(myproc == 0) then
! MS: Old version had normalization of any component equal to 
! vbin/vmax = vbin/sqrt(c^2).  Changed to make it 1.
    tot=0.; do i=-vbin,vbin; tot=tot+vix(i); enddo
    if (tot .ne. 0) then
      do i=-vbin,vbin; vix(i)=vix(i)/tot; enddo
    endif
    tot=0.; do i=-vbin,vbin; tot=tot+viy(i); enddo
    if (tot .ne. 0) then
      do i=-vbin,vbin; viy(i)=viy(i)/tot; enddo
    endif
    tot=0.; do i=-vbin,vbin; tot=tot+viz(i); enddo
    if (tot .ne. 0) then
      do i=-vbin,vbin; viz(i)=viz(i)/tot; enddo
    endif  
    tot=0.; do i=0,vbin; tot=tot+ienergy(i); enddo
    if (tot .ne. 0) then
      do i=0,vbin; ienergy(i)=ienergy(i)/tot; enddo
    endif

    tot=0.; do i=-vbin,vbin; tot=tot+vex(i); enddo
    if (tot .ne. 0) then
      do i=-vbin,vbin; vex(i)=vex(i)/tot; enddo
    endif
    tot=0.; do i=-vbin,vbin; tot=tot+vey(i); enddo
    if (tot .ne. 0) then
      do i=-vbin,vbin; vey(i)=vey(i)/tot; enddo
    endif
    tot=0.; do i=-vbin,vbin; tot=tot+vez(i); enddo
    if (tot .ne. 0) then
      do i=-vbin,vbin; vez(i)=vez(i)/tot; enddo
    endif  
    tot=0.; do i=0,vbin; tot=tot+eenergy(i); enddo
    if (tot .ne. 0) then
      do i=0,vbin; eenergy(i)=eenergy(i)/tot; enddo
    endif

    write(6,*)
    write(6,"(A)") '** Ion velocity distribution: v, vx(v), vy(v), vz(v) **'
    do i=-vbin,vbin
      write(6,"(G13.5,E13.5,E13.5,E13.5)") (i*vmaxi)/vbin,vix(i),viy(i),viz(i)
    enddo
    write(6,*)
    write(6,"(A)") '** Ion energy: Kinetic energy, KE(v) **'
    do i=0,vbin
      write(6,"(G13.5,E13.5)") (i*(vmaxi/sqrt(c_2)))/vbin,ienergy(i)
    enddo

#ifndef hybrid
    write(6,*)
    write(6,"(A)") '** Electron velocity distribution: v, vx(v), vy(v), vz(v) **'
    do i=-vbin,vbin
      write(6,"(G13.5,E13.5,E13.5,E13.5)") (i*vmaxe)/vbin,vex(i),vey(i),vez(i)
    enddo
    write(6,*)
    write(6,"(A)") '** Electron energy: Kinetic energy, KE(v) **'
    do i=0,vbin
      write(6,"(G13.5,E13.5)") (i*(vmaxe/sqrt(c_2)))/vbin,eenergy(i)
    enddo
#endif
  endif

  call MPI_Finalize(mpi_err)
End Program vdist_par

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

Subroutine diagparticles(vx,vy,vz,energy,vmax,inchannel)
  Use channels
  Implicit None
  Real(kind=drk), Dimension(-vbin:vbin) :: vx,vy,vz
  Real(kind=prk), Dimension(6,maxparticles), Save :: rv
  Real(kind=prk), Dimension(6,bufsize) :: buffer
  Real(kind=drk) :: vmax,tot
  Integer :: np,pos,n_pes_in,ch,i,j,bin

  Real(kind=drk), Dimension(0:vbin) :: energy
  Integer :: bin_en
  Real(kind=drk) :: kin_en,vbyc_sq

  Integer :: inchannel,fudge

! delete velocity bins
  do i=-vbin,vbin
    vx(i)=0.; vy(i)=0.; vz(i)=0.
  enddo
  do i = 0,vbin
    energy(i) = 0.
  enddo

! read number of original PEs
  read(unit(inchannel)) n_pes_in

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
! flag in the call to vdist_par, but who wants that?
          bin_en = int(kin_en/(vmax/sqrt(c_2))*vbin+0.5)
          if (bin_en .le. vbin) energy(bin_en)=energy(bin_en)+1.
        endif
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

End
