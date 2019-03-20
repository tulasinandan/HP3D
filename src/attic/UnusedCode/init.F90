!******************************************************************************
!                3D particle code: initialization - main entry
!                            Andreas Zeiler, 1999
!
!                        LATEST CHANGE: July 16, 2001
!
!******************************************************************************

#include "param"
#define initharris 1
#define initgem 2
#define initgem2 3
#define initharris2 4
#define initelectron 5
#define initelectron2 6
#define initelectron3 7
#define inittwostream 8
#define initelectron4 9
#define initriemann 10
#define initbeam 11
#define inittest 12
#define initperturb 13
#define initperturb2 14
#define initwave 15
#define initwaveold 16
#define initwaveB 17
#define initwaveC 18
#define initwaveD 19
#define initwaveE 20
#define initwaveF 21
#define initturb 22
#define initwaveG 23
#define initwaveH 24
#define initalfvenA 25

  Use pe_env
  Implicit None
  Integer mpi_err

#if (init_scheme == initharris)
  call init_harris_all()
#elif (init_scheme == initgem)
  call init_gem_all()
#elif (init_scheme == initgem2)
  call init_gem_all()
#elif (init_scheme == initharris2)
  call init_harris_all()
#elif (init_scheme == initelectron)
  call init_electron_all()
#elif (init_scheme == initelectron2)
  call init_electron_all()
#elif (init_scheme == initelectron3)
  call init_electron_all()
#elif (init_scheme == initelectron4)
  call init_electron_all()
#elif (init_scheme == inittwostream)
  call init_twostream()
#elif (init_scheme == initriemann)
  call init_riemann()
#elif (init_scheme == initbeam)
  call init_beam()
#elif (init_scheme == inittest)
  call init_test()
#elif (init_scheme == initperturb)
  call init_perturb_all()
#elif (init_scheme == initperturb2)
  call init_perturb_all()
#elif (init_scheme == initwave)
  call init_wave()
#elif (init_scheme == initwaveold)
  call init_waveold()
#elif (init_scheme == initwaveB)
  call init_waveB()
#elif (init_scheme == initwaveC)
  call init_waveC()
#elif (init_scheme == initwaveD)
  call init_waveD()
#elif (init_scheme == initwaveE)
  call init_waveE()
#elif (init_scheme == initwaveF)
  call init_waveF()
#elif (init_scheme == initturb)
  call init_turb()
#elif (init_scheme == initwaveG)
  call init_waveG()
#elif (init_scheme == initwaveH)
  call init_waveH()
#elif (init_scheme == initalfvenA)
  call init_alfvenA()
#else
  we never should get here
#endif
! print load balance statistics
  call load_balance()
  call MPI_Finalize(mpi_err)
  end
