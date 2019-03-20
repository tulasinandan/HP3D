!******************************************************************************
!          3D particle code: initialization - main entry, hybrid version
!                            Andreas Zeiler, 1999
!
!                         LATEST CHANGE: May 3, 2001
!
!******************************************************************************

#include "param"
#define initharris 1
#define initalfven 2
#define initenertia 3
#define initOTV 4
#define initnone 5
#define initTurb 6
#define initTwoStream 7
#define initsuperAlfvenic 8
#define initHK 9
#define initanydist 10
#define initTGV 11
#define initOTV3D 12

  Use pe_env
  Implicit None
  Integer mpi_err

#if (init_scheme == initharris-h)
  call init_harris_hyb()
#elif (init_scheme == initalfven-h)
  call init_alfven_hyb()
#elif (init_scheme == initenertia-h)
  call init_enertia_hyb()
#elif (init_scheme == initOTV-h)
  call init_OTV_hyb()
#elif (init_scheme == initnone)
  call init_none()
#elif (init_scheme == initTurb)
  call init_Turb()
#elif (init_scheme == initTwoStream)
  call init_TwoStream()
#elif (init_scheme == initsuperAlfvenic)
  call init_superAlfvenic()
#elif (init_scheme == initHK)
  call init_HK()
#elif (init_scheme == initanydist)
  call init_anydist
#elif (init_scheme == initTGV)
  call init_TGV_hyb()
#elif (init_scheme == initOTV3D)
  call init_OTV3D()
#else
  we never should get here
#endif
! print load balance statistics
  call load_balance()
  call MPI_Finalize(mpi_err)
  end
