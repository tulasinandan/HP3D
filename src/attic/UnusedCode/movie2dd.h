!******************************************************************************
!             3D particle code: header file for 2-D movie output
!                           Andreas Zeiler, 2000
!
!******************************************************************************
! M. Shay: 3/19/2006:       Added distribution function outputs.
! M. Shay: 7/14 or so/2004: Added output of ion currents.
! M. Shay: 7/28/2004:       Fixed error: ion currents where opposite
!                           sign of what they were supposed to be. 
! M. Shay: 9/14/2004:       Added output of pressure.

#ifdef movie_variables
! define variables that are needed for the movie generation
  Real(kind=drk) :: t_movieout
#ifdef double_byte
  Integer, Parameter :: mrecl = 2*nx
#else
  Integer, Parameter :: mrecl = nx
#endif
  Integer :: fileinfo,filetype,i_movie,moviecount &
            ,fileinfo2i,filetype2i,movie2i,moviecount2i  &
            ,fileinfo2e,filetype2e,movie2e,moviecount2e 
  Integer, Dimension(3) :: sizes,subsizes,starts
  Integer, Dimension(2) :: sizes2i,subsizes2i,starts2i &
         , sizes2e, subsizes2e, starts2e
  Integer, Dimension(24) :: fh
  Integer :: fh2i, fh2e
  Integer(kind=MPI_OFFSET_KIND) :: disp,offset,disp2i,disp2e & 
        ,offset2i,offset2e
#endif

#ifdef movie_open
! open movie files

  if(myproc == 0) then
    open(unit=10,file='log',action='write',form='formatted')
    open(unit=50,file='logdist',action='write',form='formatted')
  endif

! Initialize movie counter
  moviecount=0
  moviecount2i=0
  moviecount2e=0

! Define three arrays so that MPI can create the correct movie template.
! Normal movie files (fields)
  sizes = (/ pex*nx,pey*ny,pez*nz /)
  subsizes=(/ nx,ny,nz /)
  starts=(/my_pex*nx,my_pey*ny,my_pez*nz/)

! Ion distribution function
  sizes2i = (/ 2*nx*pex,nvxi /)
  subsizes2i = (/ 2*nx,nvxi /)
  starts2i = (/ 2*my_pex*nx,0 /)

! Electron distribution function
  sizes2e = (/ 2*nx*pex,nvxe /)
  subsizes2e = (/ 2*nx,nvxe /)
  starts2e = (/ 2*my_pex*nx,0 /)

! File type for usual 
  call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
				MPI_CHARACTER,filetype,mpi_err)
  call MPI_TYPE_COMMIT(filetype,mpi_err)
  call MPI_INFO_CREATE(fileinfo,mpi_err)

  call MPI_TYPE_CREATE_SUBARRAY(2,sizes2i,subsizes2i,starts2i,MPI_ORDER_FORTRAN, &
				MPI_CHARACTER,filetype2i,mpi_err)
  call MPI_TYPE_COMMIT(filetype2i,mpi_err)
  call MPI_INFO_CREATE(fileinfo2i,mpi_err)

  call MPI_TYPE_CREATE_SUBARRAY(2,sizes2e,subsizes2e,starts2e,MPI_ORDER_FORTRAN, &
				MPI_CHARACTER,filetype2e,mpi_err)
  call MPI_TYPE_COMMIT(filetype2e,mpi_err)
  call MPI_INFO_CREATE(fileinfo2e,mpi_err)



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
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'ni',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(18),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jix',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(19),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jiy',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(20),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'jiz',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(21),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'pixx',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(22),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'piyy',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(23),mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'pizz',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo,fh(24),mpi_err)

  call MPI_FILE_OPEN(MPI_COMM_WORLD,'disti',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo2i,fh2i,mpi_err)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'diste',IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
		     fileinfo2i,fh2e,mpi_err)


! Added at the recommendation of NERSCs Jonathan Carter
  disp = sizes(1)*sizes(2)*sizes(3)
  disp2i=sizes2i(1)*sizes2i(2)
  disp2e=sizes2e(1)*sizes2e(2)
  do i_movie = 1,24
    call MPI_FILE_PREALLOCATE(fh(i_movie), disp, mpi_err)
  end do
    call MPI_FILE_PREALLOCATE(fh2i, disp2i, mpi_err)
    call MPI_FILE_PREALLOCATE(fh2e, disp2e, mpi_err)
  
! The displacement must be of the type specified in the header.
  disp = 0 ; disp2i=0 ; disp2e = 0
! Apply the template to each movie file.
  do i_movie = 1,24
    call MPI_FILE_SET_VIEW(fh(i_movie),disp,MPI_CHARACTER,filetype,&
			   'native',fileinfo,mpi_err)
  enddo
  call MPI_FILE_SET_VIEW(fh2i,disp2i,MPI_CHARACTER,filetype2i,&
			   'native',fileinfo2i,mpi_err)
  call MPI_FILE_SET_VIEW(fh2e,disp2e,MPI_CHARACTER,filetype2e,&
			   'native',fileinfo2e,mpi_err)


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
    offset2i = 2*nx*nvxi*moviecount2i
    offset2e = 2*nx*nvxe*moviecount2e
    moviecount=moviecount+1
    moviecount2i=moviecount2i+1
    moviecount2e=moviecount2e+1

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
!   Output electron info
    call calc_rho(rve,np_e,real(-1.,drk),.true.)
    call calc_j(rve,np_e,real(-1.,drk),.true.)
    call calc_p(rve,np_e,real(m_e,drk),real(-1.,drk),rho,jx,jy,jz,.true.)
!print *,'myproc = ',myproc,'***',shape(pxx),'***',pxx(:,:,1)
    call volume_z(-1.*rho,fh(11),10,offset)
    call volume_z(jx,fh(12),10,offset)
    call volume_z(jy,fh(13),10,offset)
    call volume_z(jz,fh(14),10,offset)
    call volume_z(pxx,fh(15),10,offset)
    call volume_z(pyy,fh(16),10,offset)
    call volume_z(pzz,fh(17),10,offset)
!   Output ion info
    call calc_rho(rvi,np_i,real(+1.,drk),.true.)
    call calc_j(rvi,np_i,real(+1.,drk),.true.)
    call calc_p(rvi,np_i,real(1.,drk),real(1.,drk),rho,jx,jy,jz,.true.)


    call volume_z(rho,fh(18),10,offset)
    call volume_z(jx,fh(19),10,offset)
    call volume_z(jy,fh(20),10,offset)
    call volume_z(jz,fh(21),10,offset)
    call volume_z(pxx,fh(22),10,offset)
    call volume_z(pyy,fh(23),10,offset)
    call volume_z(pzz,fh(24),10,offset)
!   Reload current for stepping forward fields
    call calc_j(rvi,np_i,real(+1.,drk),.true.)
    call calc_j(rve,np_e,real(-1.,drk),.false.)
    call smooth_j(jx,jy,jz,jx,jy,jz)

    call output_dist(disti,nvxi,fh2i,50,offset2i)
    call output_dist(diste,nvxe,fh2e,50,offset2e)

! An XLF routine that flushes the I/O buffer. This way on a code crash
! the log file still has something useful in it.
    call flush_(10)
    call flush_(50)

    t_movieout=t+sign(movieout,dt)
  endif
#endif

#ifdef movie_close
! close movie files
  if (myproc == 0) close(10)
  if (myproc == 0) close(50)

  do i_movie = 1,24
    call MPI_FILE_CLOSE(fh(i_movie),mpi_err)
  enddo
  call MPI_FILE_CLOSE(fh2i,mpi_err)
  call MPI_FILE_CLOSE(fh2e,mpi_err)

#endif


