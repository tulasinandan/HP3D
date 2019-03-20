!******************************************************************************
!       3D particle code: initialization - main entry, two-fluid version
!                            Andreas Zeiler, 2001
!
!                         LATEST CHANGE: May 11, 2001
!
!******************************************************************************

#include "param"
#define initharris 1

  Use pe_env
  Implicit None
  Integer mpi_err

#if (init_scheme == initharris)
  call init_harris()
#else
  we never should get here
#endif
  call MPI_Finalize(mpi_err)
  end
