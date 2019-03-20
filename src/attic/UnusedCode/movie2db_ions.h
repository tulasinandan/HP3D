!******************************************************************************
!             3D particle code: header file for 2-D movie output
!                           Andreas Zeiler, 2000
!
!                        LATEST CHANGE: May 3, 2001
!
!******************************************************************************
! M. Shay: 7/14 or so/2004: Added output of ion currents.
! M. Shay: 7/28/2004:       Fixed error: ion currents where opposite
!    sign of what they were supposed to be. 


#ifdef movie_variables
! define variables that are needed for the movie generation
  Real(kind=drk) :: t_movieout
#ifdef double_byte
  Integer, Parameter :: mrecl = 2*nx
#else
  Integer, Parameter :: mrecl = nx
#endif
  Integer :: fileinfo,filetype,i_movie,moviecount
  Integer, Dimension(3) :: sizes,subsizes,starts
  Integer, Dimension(18) :: fh
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

  call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
				MPI_CHARACTER,filetype,mpi_err)
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
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'ni',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(11),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'ne',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(12),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jex',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(13),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jey',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(14),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jez',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(15),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jix',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(16),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jiy',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(17),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jiz',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(18),mpi_err)

! Added at the recommendation of NERSCs Jonathan Carter
  disp = sizes(1)*sizes(2)*sizes(3)
  do i_movie = 1,18
    call MPI_FILE_PREALLOCATE(fh(i_movie), disp, mpi_err)
  end do

! The displacement must be of the type specified in the header.
  disp = 0
! Apply the template to each movie file.
  do i_movie = 1,18
    call MPI_FILE_SET_VIEW(fh(i_movie),disp,MPI_CHARACTER,filetype,&
			   'native',fileinfo,mpi_err)
  enddo

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

! offset advances the file pointer for each timestep.
    offset = nx*ny*nz*moviecount
    moviecount=moviecount+1

    if (myproc==0) write(6,*) 'movie output, t=',t
    call volume_z(rho,fh(1),10,offset)
    call volume_z(jx,fh(2),10,offset)
    call volume_z(jy,fh(3),10,offset)
    call volume_z(jz,fh(4),10,offset)
    call volume_z(b1x,fh(5),10,offset) 
    call volume_z(b1y,fh(6),10,offset)
!    call poisson(temp_phi,rho,real(eps,drk))
!    call bound_phi(0,0,temp_phi)
!    call volume_z(temp_phi,fh(7),10,offset)
    call volume_z(b1z,fh(7),10,offset)
    call smooth_e(e1x,e1y,e1z,esx,esy,esz)
    call volume_z(esx,fh(8),10,offset)
    call volume_z(esy,fh(9),10,offset)
    call volume_z(esz,fh(10),10,offset)
    call calc_rho(rvi,np_i,real(+1.,drk),.true.)
    call volume_z(rho,fh(11),10,offset)
    call calc_rho(rve,np_e,real(+1.,drk),.true.)
    call volume_z(rho,fh(12),10,offset)
    call calc_j(rve,np_e,real(-1.,drk),.true.)
    call volume_z(jx,fh(13),10,offset)
    call volume_z(jy,fh(14),10,offset)
    call volume_z(jz,fh(15),10,offset)
    call calc_j(rvi,np_i,real(+1.,drk),.true.)
    call volume_z(jx,fh(16),10,offset)
    call volume_z(jy,fh(17),10,offset)
    call volume_z(jz,fh(18),10,offset)
    call calc_j(rvi,np_i,real(+1.,drk),.true.)
    call calc_j(rve,np_e,real(-1.,drk),.false.)
    call smooth_j(jx,jy,jz,jx,jy,jz)

! An XLF routine that flushes the I/O buffer. This way on a code crash
! the log file still has something useful in it.
    call flush_(10)

    t_movieout=t+sign(movieout,dt)
  endif
#endif

#ifdef movie_close
! close movie files
  if (myproc == 0) close(10)

  do i_movie = 1,15
    call MPI_FILE_CLOSE(fh(i_movie),mpi_err)
  enddo
#endif
