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
  Real(kind=drk) :: t_movieout
  Integer :: recordlength=nx*ny*nz, recnum, fileunit
  Character*24 :: filename
#endif

#ifdef movie_open
! open movie files

  if(myproc == 0) then
    open(unit=10,file='BData/log',action='write',form='formatted')
  endif


  recnum = 1   ! Initialize record counter
  fileunit = 11+myproc ! Initialize unit 


  write(filename,"(a,'.',i6.6)") 'BData/movie',myproc
  open (unit=fileunit,file=filename,status='unknown',form='unformatted' &
       ,access='direct',recl=recordlength)

! initialize movie timer
  t_movieout=0.
  do while (abs(t) > abs(t_movieout))
     t_movieout = t_movieout+sign(movieout,dt)
  enddo
! Reset record counters
  call reset_rec_ctr
#endif

#ifdef movie_write
! write to movie files
  if (abs(t) >= abs(t_movieout)) then

    if (myproc==0) write(6,*) 'movie output, t=',t
    call moviesimple(rho  ,fileunit, 10, recnum)
    call moviesimple(jx   ,fileunit, 10, recnum)
    call moviesimple(jy   ,fileunit, 10, recnum)
    call moviesimple(jz   ,fileunit, 10, recnum)
    call moviesimple(b1x  ,fileunit, 10, recnum) 
    call moviesimple(b1y  ,fileunit, 10, recnum)
    call moviesimple(b1z  ,fileunit, 10, recnum)
    call moviesimple(pe1  ,fileunit, 10, recnum)
    call moviesimple(jex  ,fileunit, 10, recnum)
    call moviesimple(jey  ,fileunit, 10, recnum)
    call moviesimple(jez  ,fileunit, 10, recnum)
    call moviesimple(jtotx,fileunit, 10, recnum)
    call moviesimple(jtoty,fileunit, 10, recnum)
    call moviesimple(jtotz,fileunit, 10, recnum)
    call moviesimple(pxx  ,fileunit, 10, recnum)
    call moviesimple(pyy  ,fileunit, 10, recnum)
    call moviesimple(pzz  ,fileunit, 10, recnum)
    call moviesimple(pxy  ,fileunit, 10, recnum)
    call moviesimple(pxz  ,fileunit, 10, recnum)
    call moviesimple(pyz  ,fileunit, 10, recnum)

! An XLF routine that flushes the I/O buffer. This way on a code crash
! the log file still has something useful in it.
#   ifdef AIX
      call flush_(10)
#   endif
#   ifdef Linux
      call flush(10)
#   endif

    t_movieout=t+sign(movieout,dt)
  endif
#endif

#ifdef movie_close
! close movie files
  if (myproc == 0) close(10)

  close(fileunit)

#endif

#ifdef moviecombine_variables
! Initialize variables to describe variable info for post-processing
Integer, parameter :: nvalues=20
Character(len=10), dimension(nvalues) :: valoutput
#endif


#ifdef movie_dataout_array
! Variable info for program to post-process movie data
  valoutput = (/ 'n ','jx  ','jy  ','jz  ','bx  ','by  ','bz  ','pe  ','jex ','jey ','jez '&
                ,'jtotx ','jtoty ','jtotz ','pxx','pyy','pzz','pxz','pyz','pxy' /)
#endif
