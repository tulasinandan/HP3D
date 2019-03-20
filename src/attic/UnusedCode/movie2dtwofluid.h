!******************************************************************************
!             3D particle code: header file for 2-D movie output
!                           Andreas Zeiler, 2000
!
!                        LATEST CHANGE: May 14, 2001
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
    open(unit=11,file='n',action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=12,file='jix', action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=13,file='jiy', action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=14,file='jiz', action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=15,file='bx', action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=16,file='by', action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=17,file='bz', action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=18,file='pe', action='write',form='unformatted',access='direct',recl=mrecl)
    open(unit=19,file='pi', action='write',form='unformatted',access='direct',recl=mrecl)
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
    call moviexy(n1,1,11,10)
    call moviexy(j1x,1,12,10)
    call moviexy(j1y,1,13,10)
    call moviexy(j1z,1,14,10)
    call moviexy(b1x,1,15,10) 
    call moviexy(b1y,1,16,10)
    call moviexy(b1z,1,17,10)
    call moviexy(pe1,1,18,10)
    call moviexy(pi1,1,19,10)
    t_movieout=t+movieout
  endif
#endif

#ifdef movie_close
! close movie files
  close(10); close(11); close(12); close(13); close(14); close(15); close(16)
  close(17); close(18); close(19)
#endif
