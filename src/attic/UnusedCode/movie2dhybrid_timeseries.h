!******************************************************************************
!             3D particle code: header file for 2-D movie output
!                       for full precision time series
!                       Tulasi Nandan Parashar, 07/2008
!
!                                                  
!
!******************************************************************************
#ifdef movie_variables
! define variables that are needed for the movie generation
!  Real(kind=drk) :: t_movieout
#ifdef double_byte
  Integer, Parameter :: mreclt = 16*nx
#else
  Integer, Parameter :: mreclt = 8*nx
#endif
#endif

#ifdef movie_open
! open movie files
  if (myproc == 0) then
    open(unit=71,file='nft'    ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
    open(unit=72,file='jixft'  ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
    open(unit=73,file='jiyft'  ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
    open(unit=74,file='jizft'  ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
    open(unit=75,file='bxft'   ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
    open(unit=76,file='byft'   ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
    open(unit=77,file='bzft'   ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
!   open(unit=78,file='peft'   ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
!   open(unit=79,file='jexft'  ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
!   open(unit=80,file='jeyft'  ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
!   open(unit=81,file='jezft'  ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
!   open(unit=82,file='jtotxft',action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
!   open(unit=83,file='jtotyft',action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
!   open(unit=84,file='jtotzft',action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
!   open(unit=85,file='pxxft'  ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
!   open(unit=86,file='pyyft'  ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
!   open(unit=87,file='pzzft'  ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
!   open(unit=88,file='pxyft'  ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
!   open(unit=89,file='pxzft'  ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
!   open(unit=90,file='pyzft'  ,action='write',form='formatted')!'unformatted',access='direct',recl=mreclt)
  endif
! Reset record counters
  call reset_rec_ctr
#endif

#ifdef movie_write
! write to movie files
    if (myproc==0) write(6,*) 'full movie output, t=',t
    call moviexyfull(rho,  1,71,.true.,.true.)
    call moviexyfull(jx,   1,72,.true.,.true.)
    call moviexyfull(jy,   1,73,.true.,.true.)
    call moviexyfull(jz,   1,74,.true.,.true.)
    call moviexyfull(b1x,  1,75,.true.,.true.) 
    call moviexyfull(b1y,  1,76,.true.,.true.)
    call moviexyfull(b1z,  1,77,.true.,.true.)
!   call moviexyfull(pe1,  1,78,.true.,.true.)
!   call moviexyfull(jex,  1,79,.true.,.true.)
!   call moviexyfull(jey,  1,80,.true.,.true.)
!   call moviexyfull(jez,  1,81,.true.,.true.)
!   call moviexyfull(jtotx,1,82,.true.,.true.)
!   call moviexyfull(jtoty,1,83,.true.,.true.)
!   call moviexyfull(jtotz,1,84,.true.,.true.)
!   call moviexyfull(pxx,  1,85,.true.,.true.)
!   call moviexyfull(pyy,  1,86,.true.,.true.)
!   call moviexyfull(pzz,  1,87,.true.,.true.)
!   call moviexyfull(pxy,  1,88,.true.,.true.)
!   call moviexyfull(pxz,  1,89,.true.,.true.)
!   call moviexyfull(pyz,  1,90,.true.,.true.)

#endif

#ifdef movie_close
! close movie files
  close(71)
  close(72)
  close(73)
  close(74)
  close(75)
  close(76)
  close(77)
! close(78)
! close(79)
! close(80)
! close(81)
! close(82)
! close(83)
! close(84)
! close(85)
! close(86)
! close(87)
! close(88)
! close(89)
! close(90)
#endif
