!******************************************************************************
!             3D particle code: header file for 2-D movie output
!                           Andreas Zeiler, 2000
!
!******************************************************************************
! M. Shay: 7/14 or so/2004: Added output of ion currents.
! M. Shay: 7/28/2004:       Fixed error: ion currents where opposite
!                           sign of what they were supposed to be. 
! M. Shay: 9/14/2004:       Added output of pressure.

! M. Shay: 8/27/2008:       Complete rewrite: Added simple output where each 
!                           processor writes to its own file. 


#ifdef movie_variables
! define variables that are needed for the movie generation
  Real(kind=drk) :: t_movieout_full!, t_dispout
!  Integer :: reclrv=8*6*maxparticles, recnumdisp, fileunitdisp
#ifdef double_precision
  Integer :: recordlengthfull=8*nx*ny*nz, recnumfull, fileunitfull
#else
  Integer :: recordlengthfull=4*nx*ny*nz, recnumfull, fileunitfull
#endif
  Character*24 :: filenamefull!, filenamedisp
#endif

#ifdef movie_open
! open movie files

  recnumfull = 1   ! Initialize record counter
  fileunitfull = 933221+myproc ! Initialize unit 


  write(filenamefull,"(a,'.',i6.6)") 'FPData/moviefull',myproc
  open (unit=fileunitfull,file=filenamefull,status='unknown',form='unformatted' &
       ,access='direct',recl=recordlengthfull)

! initialize movie timer
  t_movieout_full=0.
  do while (abs(t) > abs(t_movieout_full))
     t_movieout_full = t_movieout_full+sign(movieout_full,dt)
  enddo
! Reset record counters
  call reset_rec_ctr
#endif

#ifdef movie_write
! write to movie files
  if (abs(t) >= abs(t_movieout_full)) then

    if (myproc==0) write(6,*) 'full movie output, t=',t
    call moviesimplefull(rho  ,fileunitfull, recnumfull)
    call moviesimplefull(jx   ,fileunitfull, recnumfull)
    call moviesimplefull(jy   ,fileunitfull, recnumfull)
    call moviesimplefull(jz   ,fileunitfull, recnumfull)
    call moviesimplefull(b1x  ,fileunitfull, recnumfull) 
    call moviesimplefull(b1y  ,fileunitfull, recnumfull)
    call moviesimplefull(b1z  ,fileunitfull, recnumfull)
    call moviesimplefull(pe1  ,fileunitfull, recnumfull)
    call moviesimplefull(jtotx,fileunitfull, recnumfull)
    call moviesimplefull(jtoty,fileunitfull, recnumfull)
    call moviesimplefull(jtotz,fileunitfull, recnumfull)
    call moviesimplefull(pxx  ,fileunitfull, recnumfull)
    call moviesimplefull(pyy  ,fileunitfull, recnumfull)
    call moviesimplefull(pzz  ,fileunitfull, recnumfull)
    call moviesimplefull(pxy  ,fileunitfull, recnumfull)
    call moviesimplefull(pxz  ,fileunitfull, recnumfull)
    call moviesimplefull(pyz  ,fileunitfull, recnumfull)

    t_movieout_full=t+sign(movieout_full,dt)
  endif
!  if (abs(t) >= abs(t_dispout)) then
!
!    if (myproc==0) write(6,*) 'full movie output, t=',t
!    call moviesimplefull(rho  ,fileunitfull, recnumfull)
!    t_dispout=t+sign(dispout,dt)
!  endif
#endif

#ifdef movie_close

  close(fileunitfull)

#endif

#ifdef moviecombine_variables_full
! Initialize variables to describe variable info for post-processing
Integer, parameter :: nvalues=17
Character(len=10), dimension(nvalues) :: valoutputfull
#endif


#ifdef movie_dataout_array_full
! Variable info for program to post-process movie data
  valoutputfull = (/ 'ni ','jix   ','jiy   ','jiz   ','bx   ','by   ','bz   ','pe '&
                ,'jtotx ','jtoty ','jtotz ','pxx ','pyy ','pzz ','pxz ','pyz ','pxy ' /)
#endif
