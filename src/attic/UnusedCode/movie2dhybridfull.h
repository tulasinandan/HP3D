!******************************************************************************
!             3D particle code: header file for 2-D movie output
!                           Andreas Zeiler, 2000
!
!                        LATEST CHANGE: May 3, 2001
!
!******************************************************************************

#ifdef movie_variables
! define variables that are needed for the movie generation
  Real(kind=drk) :: t_movieout_full
#ifdef double_precision
  Integer, Parameter :: mreclf = 8*nx
#else
  Integer, Parameter :: mreclf = 4*nx
#endif
#endif

#ifdef movie_open
! open movie files
  if (myproc == 0) then
    open(unit=51,file='nf'    ,action='write',form='unformatted',access='direct',recl=mreclf)
    open(unit=52,file='jixf'  ,action='write',form='unformatted',access='direct',recl=mreclf)
    open(unit=53,file='jiyf'  ,action='write',form='unformatted',access='direct',recl=mreclf)
    open(unit=54,file='jizf'  ,action='write',form='unformatted',access='direct',recl=mreclf)
    open(unit=55,file='bxf'   ,action='write',form='unformatted',access='direct',recl=mreclf)
    open(unit=56,file='byf'   ,action='write',form='unformatted',access='direct',recl=mreclf)
    open(unit=57,file='bzf'   ,action='write',form='unformatted',access='direct',recl=mreclf)
    open(unit=58,file='pef'   ,action='write',form='unformatted',access='direct',recl=mreclf)
!   open(unit=59,file='jexf'  ,action='write',form='unformatted',access='direct',recl=mreclf)
!   open(unit=60,file='jeyf'  ,action='write',form='unformatted',access='direct',recl=mreclf)
!   open(unit=61,file='jezf'  ,action='write',form='unformatted',access='direct',recl=mreclf)
    open(unit=62,file='jtotxf',action='write',form='unformatted',access='direct',recl=mreclf)
    open(unit=63,file='jtotyf',action='write',form='unformatted',access='direct',recl=mreclf)
    open(unit=64,file='jtotzf',action='write',form='unformatted',access='direct',recl=mreclf)
    open(unit=65,file='pxxf'  ,action='write',form='unformatted',access='direct',recl=mreclf)
    open(unit=66,file='pyyf'  ,action='write',form='unformatted',access='direct',recl=mreclf)
    open(unit=67,file='pzzf'  ,action='write',form='unformatted',access='direct',recl=mreclf)
    open(unit=68,file='pxyf'  ,action='write',form='unformatted',access='direct',recl=mreclf)
    open(unit=69,file='pxzf'  ,action='write',form='unformatted',access='direct',recl=mreclf)
    open(unit=70,file='pyzf'  ,action='write',form='unformatted',access='direct',recl=mreclf)
  endif
! initialize movie timer
  t_movieout_full=0.
  do while (t > t_movieout_full)
    t_movieout_full = t_movieout_full+movieout_full
  enddo
! Reset record counters
  call reset_rec_ctr
#endif

#ifdef movie_write
! write to movie files
  if (t >= t_movieout_full) then
    if (myproc==0) write(6,*) 'full movie output, t=',t
    call moviexyfull(rho,  1,51,.false.,.false.)
    call moviexyfull(jx,   1,52,.false.,.false.)
    call moviexyfull(jy,   1,53,.false.,.false.)
    call moviexyfull(jz,   1,54,.false.,.false.)
    call moviexyfull(b1x,  1,55,.false.,.false.) 
    call moviexyfull(b1y,  1,56,.false.,.false.)
    call moviexyfull(b1z,  1,57,.false.,.false.)
    call moviexyfull(pe1,  1,58,.false.,.false.)
!   call moviexyfull(jex,  1,59,.false.,.false.)
!   call moviexyfull(jey,  1,60,.false.,.false.)
!   call moviexyfull(jez,  1,61,.false.,.false.)
    call moviexyfull(jtotx,1,62,.false.,.false.)
    call moviexyfull(jtoty,1,63,.false.,.false.)
    call moviexyfull(jtotz,1,64,.false.,.false.)
    call moviexyfull(pxx,  1,65,.false.,.false.)
    call moviexyfull(pyy,  1,66,.false.,.false.)
    call moviexyfull(pzz,  1,67,.false.,.false.)
    call moviexyfull(pxy,  1,68,.false.,.false.)
    call moviexyfull(pxz,  1,69,.false.,.false.)
    call moviexyfull(pyz,  1,70,.false.,.false.)
    t_movieout_full=t+movieout_full
  endif
#endif

#ifdef movie_close
! close movie files
  close(51)  
  close(52) 
  close(53) 
  close(54) 
  close(55) 
  close(56)
  close(57)  
  close(58)
! close(59)
! close(60)
! close(61)
  close(62)
  close(63)
  close(64)
  close(65)
  close(66)
  close(67)
  close(68)
  close(69)
  close(70)
#endif
