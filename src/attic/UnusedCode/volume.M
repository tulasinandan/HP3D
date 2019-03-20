!******************************************************************************
!                   3D particle code: diagnose volume data
!                           Andreas Zeiler, 2000
!
!                       LATEST CHANGE: July 18, 2001
!
!******************************************************************************

#include "param"

#ifndef hybrid
#ifndef twofluid
#define fullparticle
#endif
#endif

#ifndef double_byte
#define MIN 0
#define STEPS 255
#else
#define MIN (-32768)
#define STEPS 65535
#endif
#ifndef m_e
#define m_e 1.
#endif

define(PARTICLES,changequote(oq,cq)
        if (myproc == 0) then
          write(6,*) $1
          write(filename,"('volume.$2.',i3.3)") count
          open(unit=8,file=filename,action='write',form='unformatted',access='direct',recl=mrecl)
          rec_ctr = 0
        endif
        call readparticles(field,$3); call volume(field)
        if (myproc == 0) close(8) changequote())

define(FIELD,changequote(oq,cq)
        if (myproc == 0) then
          write(6,*) $1
          write(filename,"('volume.$2.',i3.3)") count
          open(unit=8,file=filename,action='write',form='unformatted',access='direct',recl=mrecl)
          rec_ctr = 0
        endif
        call readfield(field); call volume(field)
        if (myproc == 0) close(8) changequote())

undefine(`include')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Module channels
  Implicit  None
  Integer, Dimension(nchannels) :: unit
End Module channels

Module movie_io
  Integer :: rec_ctr
End Module movie_io

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Module pe_env
  Implicit	None

  include "mpif.h"

  Integer, Save :: myproc, nprocs
  Integer, Save :: mpi_drk, mpi_prk

  Integer,Save	:: my_pex, my_pey, my_pez
  Integer, Parameter :: n_pes=pexv*peyv*pezv

Contains
  Subroutine init_pe_env
    Implicit None
    Integer :: mpi_err

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

    my_pex = Mod(myproc,pexv)
    my_pey = Mod(myproc,pexv*peyv)/pexv
    my_pez = myproc/(pexv*peyv)
    return
  End Subroutine init_pe_env

  Function totalmax(data)
    Implicit None
    Real(kind=drk) :: data, maxtotal, totalmax
    Integer :: mpi_err
    call MPI_Allreduce(data,maxtotal,1,mpi_drk,MPI_MAX,MPI_COMM_WORLD,mpi_err)
    totalmax=maxtotal
  End Function totalmax

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

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Module data
  Real(kind=4), Save, Dimension(0:nxv+1,0:nyv+1,0:nzv+1) :: field
  Integer(kind=4), Dimension(0:nxv+1,0:nyv+1,0:nzv+1) :: pc
  Real(kind=drk), Save :: t,n_avg
End Module data

!------------------------------------------------------------------------------

Program volumediag
  Use channels
  Use pe_env
  Use data
  Use movie_io
  Implicit None
  Character*16 :: filename, buf
  Integer :: count,ppe,n1,n2,n3,n4,n5,ch,i,pass
  Integer :: mpi_err

#ifdef double_byte
  Integer, Parameter :: mrecl = 2*nxv
#else
  Integer, Parameter :: mrecl = nxv
#endif

  call init_pe_env()

  if (nprocs .ne. n_pes .and. myproc == 0) then
    write(6,*) '***** volume: compiled for different number of PEs *****'
    call exitallpes()
  endif

! read number of dump-up file
  call our_getarg(1,buf)
  if (buf == ' ') then
    count = 0
  else
    read (buf,*) count
  endif

  if (myproc == 0) then 
    write(filename,"(a,'-xxx.',i3.3)") head, count
    write(6,*) 'diagnose ', filename
  endif

#ifdef fullparticle
#ifndef velocorr
! read file five times
  do pass=1,5
#else
! read file ten times
  do pass=1,10
#endif
! open files
    if (myproc == 0) then
      do ch=1,nchannels
        unit(ch)=30+ch
        write(filename,"(a,'-',i3.3,'.',i3.3)") head, ch, count
        open(unit=unit(ch),file=filename,action='read',form='unformatted')
        read(unit(ch)) t,n_avg,n1,n2,n3,n4,n5
        if (n1 .ne. nxv*pexv .or. n2 .ne. nyv*peyv .or. n3 .ne. nzv*pezv &
                           .or. n4 .ne. bufsize .or. n5 .ne. nchannels) then
          write(6,*) '***** input: input-file incompatible *****'
          call exitallpes()
        endif
      enddo
    endif
    call MPI_Bcast(t,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
    call MPI_Bcast(n_avg,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
    if (pass == 1 ) then
! read electric and magnetic field and write volume data
      FIELD("Electric field E_x:",ex)
      FIELD("Electric field E_y:",ey)
      FIELD("Electric field E_z:",ez)
      FIELD("Magnetic field B_x:",bx)
      FIELD("Magnetic field B_y:",by)
      FIELD("Magnetic field B_z:",bz)
! write particle density
      PARTICLES("Ion density:",ni,0)
      PARTICLES("Electron density:",ne,0)
    else
! skip electric and magnetic field
      call skipfield()
      call skipfield()
      call skipfield()
      call skipfield()
      call skipfield()
      call skipfield()
      if (pass == 2) then
! write electron and ion velocity in x
        PARTICLES("Ion velocity - x direction:",vix,4)
        PARTICLES("Electron velocity - x direction:",vex,4)
      endif
      if (pass == 3) then
! write electron and ion velocity in y
        PARTICLES("Ion velocity - y direction:",viy,5)
        PARTICLES("Electron velocity - y direction:",vey,5)
      endif
      if (pass == 4) then
! write electron and ion velocity in z
        PARTICLES("Ion velocity - z direction:",viz,6)
        PARTICLES("Electron velocity - z direction:",vez,6)
      endif
#ifndef velocorr
! write electron and ion kinetic energy
      if (pass == 5) then
        PARTICLES("Kinetic energy per ion:",ekini,7)
        PARTICLES("Kinetic energy per electron:",ekine,8)
      endif
#else
! calculate velocity correlations
      if (pass == 5) then
        PARTICLES("Ion velocity correlation: <v_x*v_x>",vixx,9)
        PARTICLES("Electron velocity correlation: <v_x*v_x>",vexx,9)
      endif
      if (pass == 6) then
        PARTICLES("Ion velocity correlation: <v_y*v_y>",viyy,10)
        PARTICLES("Electron velocity correlation: <v_y*v_y>",veyy,10)
      endif
      if (pass == 7) then
        PARTICLES("Ion velocity correlation: <v_z*v_z>",vizz,11)
        PARTICLES("Electron velocity correlation: <v_z*v_z>",vezz,11)
      endif
      if (pass == 8) then
        PARTICLES("Ion velocity correlation: <v_x*v_y>",vixy,12)
        PARTICLES("Electron velocity correlation: <v_x*v_y>",vexy,12)
      endif
      if (pass == 9) then
        PARTICLES("Ion velocity correlation: <v_x*v_z>",vixz,13)
        PARTICLES("Electron velocity correlation: <v_x*v_z>",vexz,13)
      endif
      if (pass == 10) then
        PARTICLES("Ion velocity correlation: <v_y*v_z>",viyz,14)
        PARTICLES("Electron velocity correlation: <v_y*v_z>",veyz,14)
      endif
#endif
    endif
! close files
    do ch=1,nchannels
      if (myproc == 0) close(unit(ch))
    enddo
  enddo
#endif

#ifdef hybrid
! read file five times
  do pass=1,5
! open files
    if (myproc == 0) then
      do ch=1,nchannels
        unit(ch)=30+ch
        write(filename,"(a,'-',i3.3,'.',i3.3)") head, ch, count
        open(unit=unit(ch),file=filename,action='read',form='unformatted')
        read(unit(ch)) t,n_avg,n1,n2,n3,n4,n5
        if (n1 .ne. nxv*pexv .or. n2 .ne. nyv*peyv .or. n3 .ne. nzv*pezv &
                           .or. n4 .ne. bufsize .or. n5 .ne. nchannels) then
          write(6,*) '***** input: input-file incompatible *****'
          call exitallpes()
        endif
      enddo
    endif
    call MPI_Bcast(t,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
    call MPI_Bcast(n_avg,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
    if (pass == 1 ) then
! read magnetic field and electron pressure and write volume data
      FIELD("Magnetic field B_x:",bx)
      FIELD("Magnetic field B_y:",by)
      FIELD("Magnetic field B_z:",bz)
      FIELD("Electron pressure:",pe)
! write particle density
      PARTICLES("Ion density:",ni,0)
      if (myproc == 0) close(8)
    else
! skip electric and magnetic field
      call skipfield()
      call skipfield()
      call skipfield()
      call skipfield()
      if (pass == 2) then
! write ion velocity in x
        PARTICLES("Ion velocity - x direction:",vix,4)
      endif
      if (pass == 3) then
! write ion velocity in y
        PARTICLES("Ion velocity - y direction:",viy,5)
      endif
      if (pass == 4) then
! write ion velocity in z
        PARTICLES("Ion velocity - z direction:",viz,6)
      endif
      if (pass == 5) then
! write ion kinetic energy
        PARTICLES("Kinetic energy per ion:",ekini,7)
      endif
    endif
! close files
    do ch=1,nchannels
      if (myproc == 0) close(unit(ch))
    enddo
  enddo
#endif

#ifdef twofluid
! open files
  if (myproc == 0) then
    do ch=1,nchannels
      unit(ch)=30+ch
      write(filename,"(a,'-',i3.3,'.',i3.3)") head, ch, count
      open(unit=unit(ch),file=filename,action='read',form='unformatted')
      read(unit(ch)) t,n_avg,n1,n2,n3,n4,n5
      if (n1 .ne. nxv*pexv .or. n2 .ne. nyv*peyv .or. n3 .ne. nzv*pezv &
                         .or. n5 .ne. nchannels) then
        write(6,*) '***** input: input-file incompatible *****'
        call exitallpes()
      endif
    enddo
  endif
  call MPI_Bcast(t,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
  call MPI_Bcast(n_avg,1,mpi_drk,0,MPI_COMM_WORLD,mpi_err)
  FIELD("Magnetic field B_x:",bx)
  FIELD("Magnetic field B_y:",by)
  FIELD("Magnetic field B_z:",bz)
  FIELD("Electron pressure:",pe)
  FIELD("Ion pressure:",pi)
  FIELD("Density:",n)
  FIELD("Ion current - x direction:",jix)
  FIELD("Ion current - y direction:",jiy)
  FIELD("Ion current - z direction:",jiz)
#endif

End Program volumediag

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Subroutine readfield(f)
  Use channels
  Use pe_env
  Implicit None
  Integer x,y,z,px,py,pz,ch,npe
  Real(kind=4), Dimension(0:nxv+1,0:nyv+1,0:nzv+1) :: f
  Real(kind=4), Dimension(1:nxv*pexv) :: t
  Real(kind=drk), Dimension(1:nxv*pexv) :: input
  Integer :: mpi_err
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status

  do pz=0,pezv-1; do z=1,nzv
    do py=0,peyv-1; do y=1,nyv
      ch=mod(pz*nzv+z-1,nchannels)+1
      if (myproc == 0) then
        read (unit(ch)) input
        t=input
        do px=0,pexv-1
          npe = px+pexv*py+pexv*peyv*pz
          if (npe /= myproc) then
             call MPI_Ssend(t(px*nxv+1),nxv,MPI_REAL4,npe,123,MPI_COMM_WORLD,mpi_err)
          else
             f(1:nxv,y,z) = t(px*nxv+1:px*nxv+nxv)
          endif
        enddo
      elseif (my_pez==pz .and. my_pey==py) then
        call MPI_Recv(f(1,y,z),nxv,MPI_REAL4,0,123,MPI_COMM_WORLD,mpi_status,mpi_err)
      endif
    enddo; enddo
  enddo; enddo
  return
End Subroutine readfield

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Subroutine skipfield()
  Use channels
  Use pe_env
  Implicit None
  Integer x,y,z,px,py,pz,ch
  Real(kind=drk), Dimension(1:nxv*pexv) :: t

  do pz=0,pezv-1; do z=1,nzv
    do py=0,peyv-1; do y=1,nyv
      ch=mod(pz*nzv+z-1,nchannels)+1
      if (myproc == 0) read (unit(ch)) t
    enddo; enddo
  enddo; enddo
  return
End Subroutine skipfield

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Subroutine readparticles(f,flag)
! flag=0: density; flag=4/5/6: velocity in x/y/z
! flag=7: ion kinetic energy; flag=8: electron kinetic energy
! velocity correlations:
! flag= 9->vxx, 10->vyy, 11->vzz, 12->vxy, 13->vxz, 14->vyz
  Use channels
  Use pe_env
  Use data
  Implicit None
  Real(kind=4), Dimension(0:nxv+1,0:nyv+1,0:nzv+1) :: f
  Real(kind=prk), Dimension(6,bufsize) :: buffer
  Real(kind=drk), Parameter :: dx=lx/(nx*pex), dy=ly/(ny*pey), dz=lz/(nz*pez)
  Integer :: np,npp,n_pes_in
  Integer :: ppe, i, ch, x, y, z, flag
  Integer :: mpi_err

! delete buffer
  do z=1,nzv; do y=1,nyv; do x=1,nxv
    f(x,y,z)=0.
    pc(x,y,z)=0
  enddo; enddo; enddo

! read particles
  do ch=1,nchannels
    if (myproc == 0) read(unit(ch)) n_pes_in
  enddo
  call MPI_Bcast(n_pes_in,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)

  do ppe=0,n_pes_in-1
    ch=mod(ppe,nchannels)+1
    if (myproc == 0) then
      read(unit(ch)) npp
    endif
    call MPI_Bcast(npp,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
    do while (npp>0)
      np=0
      if (myproc == 0) then
        read(unit(ch)) buffer
        if (npp>bufsize) then
          np=bufsize
        else
          np=npp
        endif
      endif
      if (n_pes .gt. 1) then
        call MPI_Bcast(np,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
        call MPI_Bcast(buffer,6*np,mpi_prk,0,MPI_COMM_WORLD,mpi_err)
      endif
      npp=npp-bufsize
! evaluate particle data
      do i=1,np
        x=int(buffer(1,i)/dx)-my_pex*nxv+1
        if(x<1.or.x>nxv) cycle
        y=int(buffer(2,i)/dy)-my_pey*nyv+1
        if(y<1.or.y>nyv) cycle
        z=int(buffer(3,i)/dz)-my_pez*nzv+1
        if(z<1.or.z>nzv) cycle
        if (flag == 0) then 
          f(x,y,z)=f(x,y,z)+1/n_avg
        else if (flag == 7) then
          f(x,y,z)=f(x,y,z)+(buffer(4,i)**2+buffer(5,i)**2+buffer(6,i)**2)/2
          pc(x,y,z)=pc(x,y,z)+1
        else if (flag == 8) then
          f(x,y,z)=f(x,y,z)+ &
                   (buffer(4,i)**2+buffer(5,i)**2+buffer(6,i)**2)/2*m_e
          pc(x,y,z)=pc(x,y,z)+1
        else if (flag == 9) then
          f(x,y,z)=f(x,y,z)+buffer(4,i)*buffer(4,i); pc(x,y,z)=pc(x,y,z)+1
        else if (flag == 10) then
          f(x,y,z)=f(x,y,z)+buffer(5,i)*buffer(5,i); pc(x,y,z)=pc(x,y,z)+1
        else if (flag == 11) then
          f(x,y,z)=f(x,y,z)+buffer(6,i)*buffer(6,i); pc(x,y,z)=pc(x,y,z)+1
        else if (flag == 12) then
          f(x,y,z)=f(x,y,z)+buffer(4,i)*buffer(5,i); pc(x,y,z)=pc(x,y,z)+1
        else if (flag == 13) then
          f(x,y,z)=f(x,y,z)+buffer(4,i)*buffer(6,i); pc(x,y,z)=pc(x,y,z)+1
        else if (flag == 14) then
          f(x,y,z)=f(x,y,z)+buffer(5,i)*buffer(6,i); pc(x,y,z)=pc(x,y,z)+1
        else
          f(x,y,z)=f(x,y,z)+buffer(flag,i)
          pc(x,y,z)=pc(x,y,z)+1
        endif
      enddo
    enddo
  enddo
  if (flag > 0) then
    do z=1,nzv; do y=1,nyv; do x=1,nxv
      if (pc(x,y,z)>0) f(x,y,z)=f(x,y,z)/pc(x,y,z)
    enddo; enddo; enddo
  endif
  return 
End Subroutine readparticles

!------------------------------------------------------------------------------
!							  volume: x,y,z volume
!------------------------------------------------------------------------------

Subroutine volume(f)
! write xyz volume
  Use pe_env
  Implicit None

  Integer x,y,z,px,py,pz,pe,zz,yy
  Real(kind=4), Dimension(0:nxv+1,0:nyv+1,0:nzv+1) :: f
  Integer, Dimension(nxv,nyv) :: plane2d
  Integer, Dimension(nxv) :: tmp
  Integer :: mpi_err
  Integer, Dimension(MPI_STATUS_SIZE) :: mpi_status

  Real(kind=drk) :: lmin, lmax

! calculate global minimum and maximum
  lmin = 1e20; lmax = -1e20
  do z=1,nzv; do y=1,nyv; do x=1,nxv
    lmin=min(lmin,f(x,y,z)); lmax=max(lmax,f(x,y,z))
  enddo; enddo; enddo
  lmax=totalmax(lmax); lmin=totalmin(lmin)
  if (abs(lmin-lmax) .lt. eps) then
    lmin=lmin-eps;lmax=lmax+eps
  endif
  if (myproc == 0) write(6,"(E14.6,E14.6)") lmin,lmax
! loop over z planes
  do z=1,nzv*pezv
    pz=(z-1)/nzv; zz=mod(z-1,nzv)+1
    if (pz == my_pez) then
      do y=1,nyv; do x=1,nxv
        plane2d(x,y)=nint((f(x,y,zz)-lmin)/(lmax-lmin)*STEPS)+MIN
      enddo; enddo
    endif
    call MPI_Barrier(MPI_COMM_WORLD,mpi_err) ! Not absolutly necessary
    do y=1,nyv*peyv
      py=(y-1)/nyv; yy=mod(y-1,nyv)+1
      do px=0,pexv-1
        pe=px+py*pexv+pz*pexv*peyv
        if (myproc == 0) then
          if (pe /= 0) then
            call MPI_Recv(tmp,nxv,MPI_INTEGER,pe,44,MPI_COMM_WORLD,mpi_status,mpi_err)
          else
            tmp(:) = plane2d(:,yy)
          endif
          call writeline(tmp,nxv)
        elseif (myproc == pe) then
          call MPI_Ssend(plane2d(1,yy),nxv,MPI_INTEGER,0,44,MPI_COMM_WORLD,mpi_err)
        endif
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD,mpi_err) ! Not absolutly necessary
  enddo
  return
End Subroutine volume

!------------------------------------------------------------------------------
!                          support routine: writeline
!------------------------------------------------------------------------------

Subroutine writeline(tmp,nn)

  Use movie_io

  Implicit None

  Integer :: x,outh,outl,nn
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

  rec_ctr = rec_ctr+1
  write(8,rec=rec_ctr) ctmp

  return
End Subroutine writeline

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
