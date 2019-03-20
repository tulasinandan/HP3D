!******************************************************************************
!             3D particle code: header file for 3-D movie output
!                           Andreas Zeiler, 2000
!
!******************************************************************************
! M. Shay: 7/14 or so/2004: Added output of ion currents.
! M. Shay: 7/28/2004:       Fixed error: ion currents where opposite
!                           sign of what they were supposed to be. 
! M. Shay: 9/14/2004:       Added output of pressure.
! M. Swisdak: 4/2007:       Added off-diagonal components, removed ion
!                           currents (for space), tidied, renamed file.
! M. Swisdak: 3/2011:       Fixed double-byte output
! M. Swisdak: 2/2014:       Added mult_species compatability   
! M. Swisdak: 8/2014:       Added heatflux vector movies
!                           Required rewrite of nmovies and counter.
 
#ifdef movie_variables
! define variables that are needed for the movie generation
  #ifdef n_movieout
    Integer :: n_movie_offset
  #else
    Real(kind=drk) :: t_movieout
  #endif

! movies per species
  #ifdef heatfluxmovies
    Integer, Parameter :: speciesmovienum = 13
  #else
    Integer, Parameter :: speciesmovienum = 10
  #endif
  #ifdef mult_species
    Integer, Parameter :: nmovies = 10 + 3*speciesmovienum
  #else
    Integer, Parameter :: nmovies = 10 + 2*speciesmovienum
  #endif

    Integer :: fileinfo,filetype,i_movie,moviecount,fh_ind
    Integer, Dimension(3) :: sizes,subsizes,starts
    Integer, Dimension(nmovies) :: fh
    Integer(kind=MPI_OFFSET_KIND) :: disp,offset
#endif

#ifdef movie_open
! open movie files

  if(myproc == 0) then
    open(unit=10,file='log',action='write',form='formatted')
  endif

! Initialize movie counter
  moviecount=0

! Define three arrays so that MPI can create the correct movie template.
  sizes = (/ pex*nx,pey*ny,pez*nz /)
  subsizes=(/ nx,ny,nz /)
  starts=(/my_pex*nx,my_pey*ny,my_pez*nz/)

  #ifdef double_byte
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
                                MPI_INTEGER2,filetype,mpi_err)
  #elif defined four_byte
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
                                MPI_REAL,filetype,mpi_err)
  #else
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
                                MPI_CHARACTER,filetype,mpi_err)
  #endif

  call MPI_TYPE_COMMIT(filetype,mpi_err)
  call MPI_INFO_CREATE(fileinfo,mpi_err)


  call MPI_FILE_OPEN(MPI_COMM_WORLD,'rho',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(1),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jx',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(2),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jy',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(3),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jz',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(4),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'bx',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(5),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'by',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(6),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'bz',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(7),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'ex',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(8),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'ey',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(9),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'ez',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(10),mpi_err)

! Electron movies
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'ne',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(11),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jex',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(12),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jey',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(13),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jez',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(14),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'pexx',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(15),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'peyy',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(16),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'pezz',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(17),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'pexy',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(18),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'peyz',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(19),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'pexz',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(20),mpi_err)
  fh_ind = 20
#ifdef heatfluxmovies
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'qex',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+1),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'qey',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+2),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'qez',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+3),mpi_err)
  fh_ind = fh_ind+3
#endif

! Ion movies
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'ni',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
                     fileinfo,fh(fh_ind+1 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jix',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
                     fileinfo,fh(fh_ind+2 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jiy',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
                     fileinfo,fh(fh_ind+3 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jiz',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
                     fileinfo,fh(fh_ind+4 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'pixx',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
                     fileinfo,fh(fh_ind+5 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'piyy',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
                     fileinfo,fh(fh_ind+6 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'pizz',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
                     fileinfo,fh(fh_ind+7 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'pixy',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
                     fileinfo,fh(fh_ind+8 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'piyz',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
                     fileinfo,fh(fh_ind+9 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'pixz',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
                     fileinfo,fh(fh_ind+10),mpi_err)
  fh_ind = fh_ind+10
#ifdef heatfluxmovies
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'qix',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+1),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'qiy',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+2),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'qiz',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+3),mpi_err)
  fh_ind = fh_ind+3
#endif

#ifdef mult_species
! Population 2 movies
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'ni2',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+1 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jix2',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+2 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jiy2',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+3 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jiz2',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+4 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'pixx2',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+5 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'piyy2',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+6 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'pizz2',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+7 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'pixy2',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+8 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'piyz2',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+9 ),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'pixz2',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+10),mpi_err)
  fh_ind = fh_ind+10
#ifdef heatfluxmovies
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'qix2',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+1),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'qiy2',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+2),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'qiz2',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(fh_ind+3),mpi_err)
  fh_ind = fh_ind+3
#endif
#endif

! Added at the recommendation of NERSCs Jonathan Carter
  disp = sizes(1)*sizes(2)*sizes(3)
  do i_movie = 1,nmovies
    call MPI_FILE_PREALLOCATE(fh(i_movie), disp, mpi_err)
  end do

! The displacement must be of the type specified in the header.
  disp = 0
! Apply the template to each movie file.
  do i_movie = 1,nmovies
    #ifdef double_byte
      call MPI_FILE_SET_VIEW(fh(i_movie),disp,MPI_INTEGER2,filetype,&
                           'native',fileinfo,mpi_err)
    #elif defined four_byte
      call MPI_FILE_SET_VIEW(fh(i_movie),disp,MPI_REAL,filetype,&
                           'native',fileinfo,mpi_err)
    #else
      call MPI_FILE_SET_VIEW(fh(i_movie),disp,MPI_CHARACTER,filetype,&
                           'native',fileinfo,mpi_err)
    #endif
  enddo

! initialize movie timer
  #ifdef n_movieout
    n = 0
    n_movie_offset = mod( nint(abs(t)/abs(dt)),n_movieout )
  #else
    t_movieout=0.
    do while (abs(t) > abs(t_movieout))
       t_movieout = t_movieout+sign(movieout,dt)
    enddo
  #endif
#endif

#ifdef movie_write
! write to movie files
  #ifdef n_movieout
    if ( mod( n+n_movie_offset ,n_movieout) == 0) then
  #else
    if (abs(t) >= abs(t_movieout)) then
  #endif

! offset advances the file pointer for each timestep.
    offset = nx*ny*nz*moviecount
    moviecount=moviecount+1

    if (myproc==0) write(6,*) 'movie output, t=',t
    call calc_rho(rvi,np_i,1,1,.true.)
    call calc_rho(rve,np_e,-1,1,.false.)
#ifdef mult_species
    call calc_rho(rvi,np_i,qi2,massi2,.false.)
#endif

    call volume_z(rho,fh(1),10,offset)
    call volume_z(curden(1,:,:,:),fh(2),10,offset)
    call volume_z(curden(2,:,:,:),fh(3),10,offset)
    call volume_z(curden(3,:,:,:),fh(4),10,offset)
    call volume_z(b1x,fh(5),10,offset) 
    call volume_z(b1y,fh(6),10,offset)
    call volume_z(b1z,fh(7),10,offset)
    call volume_z(e1x,fh(8),10,offset)
    call volume_z(e1y,fh(9),10,offset)
    call volume_z(e1z,fh(10),10,offset)

! Electron data
    call calc_rho(rve,np_e,-1,1,.true.)
    call calc_j(rve,np_e,-1,1,.true.)
    call calc_p(rve,np_e,-1,1,.true.,.false.)
! The absolute value of the charge density
    call volume_z(-rho,fh(11),10,offset)
    call volume_z(curden(1,:,:,:),fh(12),10,offset)
    call volume_z(curden(2,:,:,:),fh(13),10,offset)
    call volume_z(curden(3,:,:,:),fh(14),10,offset)
    call volume_z(pressure(1,:,:,:),fh(15),10,offset)
    call volume_z(pressure(2,:,:,:),fh(16),10,offset)
    call volume_z(pressure(3,:,:,:),fh(17),10,offset)
    call volume_z(pressure(4,:,:,:),fh(18),10,offset)
    call volume_z(pressure(5,:,:,:),fh(19),10,offset)
    call volume_z(pressure(6,:,:,:),fh(20),10,offset)
    fh_ind = 20
#ifdef heatfluxmovies
    call volume_z(heatflux(1,:,:,:),fh(fh_ind+1),10,offset)
    call volume_z(heatflux(2,:,:,:),fh(fh_ind+2),10,offset)
    call volume_z(heatflux(3,:,:,:),fh(fh_ind+3),10,offset)
    fh_ind = fh_ind+3
#endif

! Ion data
    call calc_rho(rvi,np_i,1,1,.true.)
    call calc_j(rvi,np_i,1,1,.true.)
    call calc_p(rvi,np_i,1,1,.true.,.false.)
    call volume_z(rho,fh(fh_ind+1),10,offset)
    call volume_z(curden(1,:,:,:),fh(fh_ind+2),10,offset)
    call volume_z(curden(2,:,:,:),fh(fh_ind+3),10,offset)
    call volume_z(curden(3,:,:,:),fh(fh_ind+4),10,offset)
    call volume_z(pressure(1,:,:,:),fh(fh_ind+5 ),10,offset)
    call volume_z(pressure(2,:,:,:),fh(fh_ind+6 ),10,offset)
    call volume_z(pressure(3,:,:,:),fh(fh_ind+7 ),10,offset)
    call volume_z(pressure(4,:,:,:),fh(fh_ind+8 ),10,offset)
    call volume_z(pressure(5,:,:,:),fh(fh_ind+9 ),10,offset)
    call volume_z(pressure(6,:,:,:),fh(fh_ind+10),10,offset)
    fh_ind = fh_ind+10
#ifdef heatfluxmovies
    call volume_z(heatflux(1,:,:,:),fh(fh_ind+1),10,offset)
    call volume_z(heatflux(2,:,:,:),fh(fh_ind+2),10,offset)
    call volume_z(heatflux(3,:,:,:),fh(fh_ind+3),10,offset)
    fh_ind = fh_ind+3
#endif

#ifdef mult_species
! Population 2 Ion data
    call calc_rho(rvi,np_i,qi2,massi2,.true.)
    call calc_j(rvi,np_i,qi2,massi2,.true.)
    call calc_p(rvi,np_i,qi2,massi2,.true.,.true.)
    call volume_z(rho,fh(fh_ind+1),10,offset)
    call volume_z(curden(1,:,:,:),fh(fh_ind+2),10,offset)
    call volume_z(curden(2,:,:,:),fh(fh_ind+3),10,offset)
    call volume_z(curden(3,:,:,:),fh(fh_ind+4),10,offset)
    call volume_z(pressure(1,:,:,:),fh(fh_ind+5 ),10,offset) !pixx2
    call volume_z(pressure(2,:,:,:),fh(fh_ind+6 ),10,offset) !piyy2
    call volume_z(pressure(3,:,:,:),fh(fh_ind+7 ),10,offset) !pizz2
    call volume_z(pressure(4,:,:,:),fh(fh_ind+8 ),10,offset) !pixy2
    call volume_z(pressure(5,:,:,:),fh(fh_ind+9 ),10,offset) !piyz2
    call volume_z(pressure(6,:,:,:),fh(fh_ind+10),10,offset) !pixz2
    fh_ind=fh_ind+10
#ifdef heatfluxmovies
    call volume_z(heatflux(1,:,:,:),fh(fh_ind+1),10,offset)
    call volume_z(heatflux(2,:,:,:),fh(fh_ind+2),10,offset)
    call volume_z(heatflux(3,:,:,:),fh(fh_ind+3),10,offset)
    fh_ind = fh_ind+3
#endif
#endif

!   Reload current for stepping forward fields
    call calc_j(rvi,np_i,1,1,.true.)
    call calc_j(rve,np_e,-1,1,.false.)
#ifdef mult_species
    call calc_j(rvi,np_i,qi2,massi2,.false.)
#endif

#ifdef smooth_e_j_rho
    call smooth_j(curden(1,:,:,:),curden(2,:,:,:),curden(3,:,:,:), &
                  curden(1,:,:,:),curden(2,:,:,:),curden(3,:,:,:))
#endif

! An XLF routine that flushes the I/O buffer. This way on a code crash
! the log file still has something useful in it.
    call flush(10)

    #ifdef n_movieout
    #else
      t_movieout=t+sign(movieout,dt)
    #endif
  endif
#endif

#ifdef movie_close
! close movie files
  if (myproc == 0) close(10)
  do i_movie = 1,nmovies
    call MPI_FILE_CLOSE(fh(i_movie),mpi_err)
  enddo

#endif
