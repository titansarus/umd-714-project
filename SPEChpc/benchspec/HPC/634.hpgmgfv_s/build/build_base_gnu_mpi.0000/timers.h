//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#ifndef TIMER_H
#define TIMER_H

  #include<stdint.h>

  #if defined(SPEC_OPENMP) || defined(SPEC_OPENMP_TARGET)
    #include <omp.h>
    #define getTime() (omp_get_wtime())

  #elif defined(USE_MPI)
    #include <mpi.h>
    #define getTime() (MPI_Wtime())

  #else
    // user must provide a function getTime and include it in timers.c
    // if calibration is necesary, then the user must #define CALIBRATE_TIMER
    double getTime();
  #endif

#endif
