!******************************************************************************
!             3D particle code: header file for 2-D movie output
!                           Andreas Zeiler, 2000
!
!                        LATEST CHANGE: May 3, 2001
!
!******************************************************************************

#ifdef movie_variables
! define variables that are needed for the movie generation
  Real(kind=drk) :: t_movieout
#ifdef double_byte
  Integer, Parameter :: mrecl = 2*nx
#else
  Integer, Parameter :: mrecl = nx
#endif
#endif

#ifdef movie_open
! open movie files
  if (myproc == 0) then
    open(unit=10,file='log',action='write',form='formatted')
    open(unit=11,file='n',  action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=12,file='jix',action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=13,file='jiy',action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=14,file='jiz',action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=15,file='bx', action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=16,file='by', action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=17,file='bz', action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=18,file='pe', action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=27,file='jex',action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=28,file='jey',action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=29,file='jez',action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=30,file='jtotx',action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=31,file='jtoty',action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=32,file='jtotz',action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=33,file='pxx',action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=34,file='pyy',action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=35,file='pzz',action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=36,file='pxy',action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=37,file='pxz',action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=38,file='pyz',action='write',form='unformatted',access='direct',recl=mrecl)
  endif
! initialize movie timer
  t_movieout=0.
  do while (t > t_movieout)
    t_movieout = t_movieout+movieout
  enddo
! Reset record counters
  call reset_rec_ctr
#endif

#ifdef movie_write
! write to movie files
  if (t >= t_movieout) then
    if (myproc==0) write(6,*) 'movie output, t=',t
    call moviexy(rho,1,11,10)
    call moviexy(jx,1,12,10)
    call moviexy(jy,1,13,10)
    call moviexy(jz,1,14,10)
    call moviexy(b1x,1,15,10) 
    call moviexy(b1y,1,16,10)
    call moviexy(b1z,1,17,10)
    call moviexy(pe1,1,18,10)
    call moviexy(jex,1,27,10)
    call moviexy(jey,1,28,10)
    call moviexy(jez,1,29,10)
    call moviexy(jtotx,1,30,10)
    call moviexy(jtoty,1,31,10)
    call moviexy(jtotz,1,32,10)
    call moviexy(pxx,1,33,10)
    call moviexy(pyy,1,34,10)
    call moviexy(pzz,1,35,10)
    call moviexy(pxy,1,36,10)
    call moviexy(pxz,1,37,10)
    call moviexy(pyz,1,38,10)
    t_movieout=t+movieout

  endif
#endif

#ifdef movie_close
! close movie files
  close(10); close(11); close(12); close(13); close(14); close(15); close(16)
  close(17); close(18); close(27); close(28); close(29); close(30); close(31)
  close(32); close(33); close(34); close(35); close(36); close(37); close(38)
#endif
