//------------------------------------------------------------------------------------------------------------------------------
// Copyright Notice 
//------------------------------------------------------------------------------------------------------------------------------
// HPGMG, Copyright (c) 2014, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to
// receipt of any required approvals from the U.S. Dept. of Energy).  All
// rights reserved.
// 
// If you have questions about your rights to use or distribute this
// software, please contact Berkeley Lab's Technology Transfer Department
// at  TTD@lbl.gov.
// 
// NOTICE.  This software is owned by the U.S. Department of Energy.  As
// such, the U.S. Government has been granted for itself and others
// acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide
// license in the Software to reproduce, prepare derivative works, and
// perform publicly and display publicly.  Beginning five (5) years after
// the date permission to assert copyright is obtained from the U.S.
// Department of Energy, and subject to any subsequent five (5) year
// renewals, the U.S. Government is granted for itself and others acting
// on its behalf a paid-up, nonexclusive, irrevocable, worldwide license
// in the Software to reproduce, prepare derivative works, distribute
// copies to the public, perform publicly and display publicly, and to
// permit others to do so.
//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#ifdef VERBOSE
# define _GNU_SOURCE /* needed for sched_getaffinity */
# include <sched.h> /* defines cpu_set_t and sched_getaffinity */
# include <sys/types.h> /* defines pid_t */
# include <unistd.h> /* defines syscall */
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
//------------------------------------------------------------------------------------------------------------------------------
#ifdef USE_MPI
#include <mpi.h>
#endif
#if defined(SPEC_OPENMP) || defined(SPEC_OPENMP_TARGET)
#include <omp.h>
#endif
#if defined(SPEC_OPENACC)
#include <openacc.h>
#endif
#ifdef PROFILE_SOLVE_ONLY
#include <cuda_profiler_api.h>
#endif
#ifdef SPEC
#include <specmpitime.h>
#endif
//------------------------------------------------------------------------------------------------------------------------------
#include "timers.h"
#include "defines.h"
#include "level.h"
#include "mg.h"
#include "operators.h"
#include "solvers.h"
#include "cuda/common.h"
#include "offload-fns.h"

#ifdef VERBOSE
void show_environment(int mpi);
void show_process_environment(int mpi);
void show_process_affinity(int mpi);
#if defined(SPEC_OPENMP) || defined(SPEC_OPENMP_TARGET)
void show_omp_thread_affinity(int mpi);
# endif
static char *cpuset_to_cstr(cpu_set_t *mask, char *str);
void show_procmem(void);
double get_procmem(void);
#endif

//------------------------------------------------------------------------------------------------------------------------------
// Print device info
int cudaCheckPeerToPeer(int rank){
  int ndev = 0;

#ifdef SPEC_CUDA
  int peer = 0;
  // query number of GPU devices in the system
  cudaGetDeviceCount(&ndev);
  // print only for 10 first ranks
  // if (rank < 10) printf("rank %d:  Number of visible GPUs:  %d\n",rank,ndev); // Too verbose at scale

  // check for peer to peer mappings
  // print only for rank 0
  for(int i=0;i<ndev;i++)
  for(int j=i+1;j<ndev;j++){
    struct cudaDeviceProp devPropi,devPropj;
    cudaGetDeviceProperties(&devPropi,i);
    cudaGetDeviceProperties(&devPropj,j);
    cudaDeviceCanAccessPeer(&peer,i,j);
    // this info can also be collected with nvidia-smi topo -m
    //printf("rank %d:  Peer access from %s (GPU%d) -> %s (GPU%d) : %s\n",rank,devPropi.name,i,devPropj.name,j,peer?"Yes":"No"); // Too verbose at scale
  }
#endif
  return ndev;
}

void bench_hpgmg(mg_type *all_grids, int onLevel, int solves_per_level, double a, double b, double dtol, double rtol){

  #ifdef USE_HPM // IBM performance counters for BGQ...
  if( onLevel==0 )HPM_Start("FMGSolve()");
  #endif

  if(all_grids->levels[onLevel]->my_rank==0){
    fprintf(stdout,"\n\n===== Running %d solves ========================================================\n",solves_per_level);
    fflush(stdout);
  }

  int numSolves =  0; // solves completed
  MGResetTimers(all_grids);
  while( (numSolves<solves_per_level) ){
    zero_vector(all_grids->levels[onLevel],VECTOR_U);
    #ifdef USE_FCYCLES
    FMGSolve(all_grids,onLevel,VECTOR_U,VECTOR_F,a,b,dtol,rtol);
    #else
    MGSolve(all_grids,onLevel,VECTOR_U,VECTOR_F,a,b,dtol,rtol);
    #endif
    numSolves++;
  }

  #ifdef USE_HPM // IBM performance counters for BGQ...
  if( onLevel==0 )HPM_Stop("FMGSolve()");
  #endif
}


//------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[]){
  int my_rank=0;
  int num_tasks=1;
#ifdef VERBOSE
  int OMP_Threads = 1; // Only used in a print statement
# ifdef SPEC_OPENMP
# pragma omp parallel
  {
# pragma omp master
    {
      OMP_Threads = omp_get_num_threads();
    }
  }
# endif
#endif

#if defined(VERBOSE) || (defined(USE_MPI) && \
			 (defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET) || defined(SPEC_CUDA)))
  int num_devices = 0; // Used in a print statement but also to assign MPI ranks to devices
#endif

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  // initialize MPI and HPM
  #ifdef USE_MPI
  int    actual_threading_model = -1;
  int requested_threading_model = -1;
      requested_threading_model = MPI_THREAD_SINGLE;
    //requested_threading_model = MPI_THREAD_FUNNELED;
    //requested_threading_model = MPI_THREAD_SERIALIZED;
    //requested_threading_model = MPI_THREAD_MULTIPLE;
    #ifdef SPEC_OPENMP
      requested_threading_model = MPI_THREAD_FUNNELED;
    //requested_threading_model = MPI_THREAD_SERIALIZED;
    //requested_threading_model = MPI_THREAD_MULTIPLE;
    #endif
    #ifdef USE_MPI_THREAD_MULTIPLE
      requested_threading_model = MPI_THREAD_MULTIPLE;
    #endif
  MPI_Init_thread(&argc, &argv, requested_threading_model, &actual_threading_model);
  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
SPEC_TIME_START
SPEC_TIME_START_INIT

  // Set device for this rank...
#if defined(SPEC_OPENMP_TARGET) || defined(SPEC_OPENACC)
  int my_device;
  MPI_Comm shmcomm;
  int local_rank;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                      MPI_INFO_NULL, &shmcomm);
  MPI_Comm_rank(shmcomm, &local_rank);
# ifdef SPEC_OPENMP_TARGET
  num_devices = omp_get_num_devices();
  my_device = num_devices >= 1 ? local_rank % num_devices : omp_get_initial_device();
  omp_set_default_device(my_device);
# else //OpenACC
  acc_device_t my_device_type;
  my_device_type = acc_get_device_type();
  num_devices = acc_get_num_devices(my_device_type);
  my_device = local_rank % num_devices;
  acc_set_device_num(my_device, my_device_type);
# endif
#endif // set rank to device

  #ifdef USE_HPM // IBM HPM counters for BGQ...
  HPM_Init();
  #endif
  #endif // USE_MPI

#if defined(VERBOSE) && defined(PROFILE_SOLVE_ONLY)
  if (my_rank == 0) fprintf(stderr, "Profiling the solve phase only\n");
#endif

  check_consistent_build();
  set_cuda_execution_array();
  int rtn = device_runtime_init();
  if (rtn != 0) {
    printf("Error: the device runtime did not initialize correctly\n");
#ifdef USE_MPI
    MPI_Finalize();
#endif
    exit(1);
  }
#ifdef VERBOSE
  show_environment(my_rank);
  fflush(stdout);
  fflush(stderr);
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
  //CD: Comment out because I want more detail and don't want nested regions
  //NVTX_PUSH("main",1)  // start NVTX profiling

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  // parse the arguments...
  int     log2_box_dim           =  6; // 64^3
  // int     target_boxes_per_rank  =  1; // Original method
  // int64_t target_memory_per_rank = -1; // not specified
  int64_t box_dim                = -1;
  int64_t boxes_in_i             = -1;
  // int64_t target_boxes           = -1; // Original method
  int     log2_grid_dim;
  int     solves_per_level;

  if(argc==4){
    log2_box_dim=atoi(argv[1]);
    // target_boxes_per_rank=atoi(argv[2]); // CD: Original method
    log2_grid_dim=atoi(argv[2]); // CD: New method keeps problem size fixed
    solves_per_level=atoi(argv[3]); // CD: New method keeps number of solves fixed

    if (log2_box_dim > log2_grid_dim) {
      if(my_rank==0){fprintf(stderr,"log2_box_dim cannot be bigger than log2_grid_dim\n");}
      #ifdef USE_MPI
      MPI_Finalize();
      #endif
      exit(0);
    }

    if(log2_box_dim>9){
      // NOTE, in order to use 32b int's for array indexing, box volumes must be less than 2^31 doubles
      if(my_rank==0){fprintf(stderr,"log2_box_dim must be less than 10\n");}
      #ifdef USE_MPI
      MPI_Finalize();
      #endif
      exit(0);
    }

    if(log2_box_dim<4){
      if(my_rank==0){fprintf(stderr,"log2_box_dim must be at least 4\n");}
      #ifdef USE_MPI
      MPI_Finalize();
      #endif
      exit(0);
    }

    if(solves_per_level<1){
      if(my_rank==0){fprintf(stderr,"solves_per_level must be at least 1\n");}
      #ifdef USE_MPI
      MPI_Finalize();
      #endif
      exit(0);
    }

#if 0
    if(target_boxes_per_rank<1){
      if(my_rank==0){fprintf(stderr,"target_boxes_per_rank must be at least 1\n");}
      #ifdef USE_MPI
      MPI_Finalize();
      #endif
      exit(0);
    }

    #ifndef MAX_COARSE_DIM
    #define MAX_COARSE_DIM 11
    #endif
    box_dim=1<<log2_box_dim;
    target_boxes = (int64_t)target_boxes_per_rank*(int64_t)num_tasks;
    boxes_in_i = -1;
    int64_t bi;
    for(bi=1;bi<1000;bi++){ // search all possible problem sizes to find acceptable boxes_in_i
      int64_t total_boxes = bi*bi*bi;
      if(total_boxes<=target_boxes){
        int64_t coarse_grid_dim = box_dim*bi;
        while( (coarse_grid_dim%2) == 0){coarse_grid_dim=coarse_grid_dim/2;}
        if(coarse_grid_dim<=MAX_COARSE_DIM){
          boxes_in_i = bi;
        }
      }
    }
#endif

    box_dim=1<<log2_box_dim;
    // log2 input parameters guarantees boxes_in_i leaves no remainder
    boxes_in_i = 1 << (log2_grid_dim - log2_box_dim);

    if(boxes_in_i<1){
      if(my_rank==0){fprintf(stderr,"failed to find an acceptable problem size\n");}
      #ifdef USE_MPI
      MPI_Finalize();
      #endif
      exit(0);
    }
  } // argc==3

  #if 0
  else if(argc==2){ // interpret argv[1] as target_memory_per_rank
    char *ptr = argv[1];
    char *tmp;
    target_memory_per_rank = strtol(ptr,&ptr,10);
    if(target_memory_per_rank<1){
      if(my_rank==0){fprintf(stderr,"unrecognized target_memory_per_rank... '%s'\n",argv[1]);}
      #ifdef USE_MPI
      MPI_Finalize();
      #endif
      exit(0);
    }
    tmp=strstr(ptr,"TB");if(tmp){ptr=tmp+2;target_memory_per_rank *= (uint64_t)(1<<30)*(1<<10);}
    tmp=strstr(ptr,"GB");if(tmp){ptr=tmp+2;target_memory_per_rank *= (uint64_t)(1<<30);}
    tmp=strstr(ptr,"MB");if(tmp){ptr=tmp+2;target_memory_per_rank *= (uint64_t)(1<<20);}
    tmp=strstr(ptr,"tb");if(tmp){ptr=tmp+2;target_memory_per_rank *= (uint64_t)(1<<30)*(1<<10);}
    tmp=strstr(ptr,"gb");if(tmp){ptr=tmp+2;target_memory_per_rank *= (uint64_t)(1<<30);}
    tmp=strstr(ptr,"mb");if(tmp){ptr=tmp+2;target_memory_per_rank *= (uint64_t)(1<<20);}
    if( (ptr) && (*ptr != '\0') ){
      if(my_rank==0){fprintf(stderr,"unrecognized units... '%s'\n",ptr);}
      #ifdef USE_MPI
      MPI_Finalize();
      #endif
      exit(0);
    }
    // FIX, now search for an 'acceptable' box_dim and boxes_in_i constrained by target_memory_per_rank, num_tasks, and MAX_COARSE_DIM
  } // argc==2
  #endif


  else{
    if(my_rank==0){fprintf(stderr,"usage: ./hpgmg-fv  log2_box_dim  log2_grid_dim  solves_per_level\n");}
    // if(my_rank==0){fprintf(stderr,"usage: ./hpgmg-fv  [log2_box_dim]  [target_boxes_per_rank]\n");}
                 //fprintf(stderr,"       ./hpgmg-fv  [target_memory_per_rank[MB,GB,TB]]\n");}
    #ifdef USE_MPI
    MPI_Finalize();
    #endif
    exit(0);
  }



 
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if(my_rank==0){
  fprintf(stdout,"\n\n");
  fprintf(stdout,"********************************************************************************\n");
  fprintf(stdout,"***                            HPGMG-FV Benchmark                            ***\n");
  fprintf(stdout,"********************************************************************************\n");
#ifdef VERBOSE
  #ifdef USE_MPI
       if(requested_threading_model == MPI_THREAD_MULTIPLE  )fprintf(stdout,"Requested MPI_THREAD_MULTIPLE, ");
  else if(requested_threading_model == MPI_THREAD_SINGLE    )fprintf(stdout,"Requested MPI_THREAD_SINGLE, ");
  else if(requested_threading_model == MPI_THREAD_FUNNELED  )fprintf(stdout,"Requested MPI_THREAD_FUNNELED, ");
  else if(requested_threading_model == MPI_THREAD_SERIALIZED)fprintf(stdout,"Requested MPI_THREAD_SERIALIZED, ");
  else if(requested_threading_model == MPI_THREAD_MULTIPLE  )fprintf(stdout,"Requested MPI_THREAD_MULTIPLE, ");
  else                                                       fprintf(stdout,"Requested Unknown MPI Threading Model (%d), ",requested_threading_model);
       if(actual_threading_model    == MPI_THREAD_MULTIPLE  )fprintf(stdout,"got MPI_THREAD_MULTIPLE\n");
  else if(actual_threading_model    == MPI_THREAD_SINGLE    )fprintf(stdout,"got MPI_THREAD_SINGLE\n");
  else if(actual_threading_model    == MPI_THREAD_FUNNELED  )fprintf(stdout,"got MPI_THREAD_FUNNELED\n");
  else if(actual_threading_model    == MPI_THREAD_SERIALIZED)fprintf(stdout,"got MPI_THREAD_SERIALIZED\n");
  else if(actual_threading_model    == MPI_THREAD_MULTIPLE  )fprintf(stdout,"got MPI_THREAD_MULTIPLE\n");
  else                                                       fprintf(stdout,"got Unknown MPI Threading Model (%d)\n",actual_threading_model);
  #endif
  fprintf(stdout,"%d MPI Tasks of %d threads\n",num_tasks,OMP_Threads);
#endif
  fprintf(stdout,"\n\n===== Benchmark setup ==========================================================\n");
  }

  NVTX_PUSH("InitFineLevel",1)
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // create the fine level...
  #ifdef USE_PERIODIC_BC
  int bc = BC_PERIODIC;
  int minCoarseDim = 2; // avoid problems with black box calculation of D^{-1} for poisson with periodic BC's on a 1^3 grid
  #else
  int bc = BC_DIRICHLET;
  int minCoarseDim = 1; // assumes you can drop order on the boundaries
  #endif
  level_type level_h;
  int ghosts=stencil_get_radius();
  create_level(&level_h,boxes_in_i,box_dim,ghosts,VECTORS_RESERVED,bc,my_rank,num_tasks,NULL);
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #ifdef USE_HELMHOLTZ
  double a=1.0;double b=1.0; // Helmholtz
  if(my_rank==0)fprintf(stdout,"  Creating Helmholtz (a=%f, b=%f) test problem\n",a,b);
  #else
  double a=0.0;double b=1.0; // Poisson
  if(my_rank==0)fprintf(stdout,"  Creating Poisson (a=%f, b=%f) test problem\n",a,b);
  #endif
  double h=1.0/( (double)boxes_in_i*(double)box_dim );  // [0,1]^3 problem
  initialize_problem(&level_h,h,a,b);                   // initialize VECTOR_ALPHA, VECTOR_BETA*, and VECTOR_F

#if defined(SPEC_OPENACC) && !defined(SPEC_MANAGED_MEMORY)
  create_level_device(&level_h);
  build_exchange_ghosts_device(&level_h);
  build_box_device(&level_h,0);
#endif

#if defined(SPEC_OPENMP_TARGET) && !defined(SPEC_MANAGED_MEMORY)
  create_level_device_openmp(&level_h);
  build_exchange_ghosts_device_openmp(&level_h);
  build_box_device_openmp(&level_h,0);
#endif

  rebuild_operator(&level_h,NULL,a,b);                  // calculate Dinv and lambda_max
  if(level_h.boundary_condition.type == BC_PERIODIC){   // remove any constants from the RHS for periodic problems
    double average_value_of_f = mean(&level_h,VECTOR_F);
    if(average_value_of_f!=0.0){
      if(my_rank==0){fprintf(stderr,"  WARNING... Periodic boundary conditions, but f does not sum to zero... mean(f)=%e\n",average_value_of_f);}
      shift_vector(&level_h,VECTOR_F,VECTOR_F,-average_value_of_f);
    }
  }
  // End of NVTX region: InitFineLevel
  NVTX_POP
SPEC_TIME_STOP_INIT
SPEC_TIME_START_COMP
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // create the MG hierarchy...
  mg_type MG_h;
  NVTX_PUSH("MGBuild",2)
  MGBuild(&MG_h,&level_h,a,b,minCoarseDim);             // build the Multigrid Hierarchy
  NVTX_POP

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // HPGMG-500 benchmark proper
  // evaluate performance on problem sizes of h, 2h, and 4h
  // (i.e. examine dynamic range for problem sizes N, N/8, and N/64)
//double dtol=1e-15;double rtol=  0.0; // converged if ||D^{-1}(b-Ax)|| < dtol
  double dtol=  0.0;double rtol=1e-10; // converged if ||b-Ax|| / ||b|| < rtol
  int l;
  #ifndef TEST_ERROR
  #ifdef VERBOSE
  double AverageSolveTime[3];
  #endif
  #ifdef USE_NVTX
  const char* nvtx_solve_id_label[3] = {"MGSolve-0", "MGSolve-1", "MGSolve-2"};
  #endif
  #ifdef PROFILE_SOLVE_ONLY
  cudaProfilerStart();
  #endif
  for(l=0;l<3;l++){
    NVTX_PUSH(nvtx_solve_id_label[l],3+l)
    if(l>0)restriction(MG_h.levels[l],VECTOR_F,MG_h.levels[l-1],VECTOR_F,RESTRICT_CELL);
    bench_hpgmg(&MG_h,l,solves_per_level,a,b,dtol,rtol);
    NVTX_POP
#ifdef VERBOSE
    AverageSolveTime[l] = (double)MG_h.timers.MGSolve / (double)MG_h.MGSolves_performed;
    if(my_rank==0){fprintf(stdout,"\n\n===== Timing Breakdown =========================================================\n");}
    MGPrintTiming(&MG_h,l);
#endif
  }
  #ifdef PROFILE_SOLVE_ONLY
  cudaProfilerStop();
  #endif
SPEC_TIME_STOP_COMP

#ifdef VERBOSE
  if(my_rank==0){
    #ifdef CALIBRATE_TIMER
    double _timeStart=getTime();sleep(1);double _timeEnd=getTime();
    double SecondsPerCycle = (double)1.0/(double)(_timeEnd-_timeStart);
    #else
    double SecondsPerCycle = 1.0;
    #endif
    fprintf(stdout,"\n\n===== Performance Summary ======================================================\n");
    for(l=0;l<3;l++){
      double DOF = (double)MG_h.levels[l]->dim.i*(double)MG_h.levels[l]->dim.j*(double)MG_h.levels[l]->dim.k;
      double seconds = SecondsPerCycle*(double)AverageSolveTime[l];
      double DOFs = DOF / seconds;
      int used_devices = (num_devices > 0 && num_tasks < num_devices) ? num_tasks : num_devices;
      fprintf(stdout,"  h=%0.15e  DOF=%0.15e  time=%0.6f  DOF/s=%0.3e  MPI=%d  OMP=%d  ACC=%d\n",
	      MG_h.levels[l]->h,DOF,seconds,DOFs,num_tasks,OMP_Threads,used_devices);
    }
  }
  #endif
#endif

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if(my_rank==0){fprintf(stdout,"\n\n===== Richardson error analysis ================================================\n");}
  // solve A^h u^h = f^h
  // solve A^2h u^2h = f^2h
  // solve A^4h u^4h = f^4h
  // error analysis...
  MGResetTimers(&MG_h);
  for(l=0;l<3;l++){
    if(l>0)restriction(MG_h.levels[l],VECTOR_F,MG_h.levels[l-1],VECTOR_F,RESTRICT_CELL);
    zero_vector(MG_h.levels[l],VECTOR_U);
    #ifdef USE_FCYCLES
    FMGSolve(&MG_h,l,VECTOR_U,VECTOR_F,a,b,dtol,rtol);
    #else
     MGSolve(&MG_h,l,VECTOR_U,VECTOR_F,a,b,dtol,rtol);
    #endif
  }
  //CD: Comment out because I want more detail and don't want nested regions
  //NVTX_POP  // stop NVTX profiling
  richardson_error(&MG_h,0,VECTOR_U);
#ifdef VERBOSE
  show_procmem();
#endif

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if(my_rank==0){fprintf(stdout,"\n\n===== Deallocating memory ======================================================\n");}
  MGDestroy(&MG_h);
  destroy_level(&level_h);


  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if(my_rank==0){fprintf(stdout,"\n\n===== Done =====================================================================\n");}
SPEC_TIME_STOP
SPEC_TIME_FINAL(true,0.0)
  #ifdef USE_MPI
  #ifdef USE_HPM // IBM performance counters for BGQ...
  HPM_Print();
  #endif
  MPI_Finalize();
  #endif
  return 0;
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
}

#ifdef VERBOSE
void show_environment(int mpi)
{
  show_process_environment(mpi);
#if defined(SPEC_OPENMP) || defined(SPEC_OPENMP_TARGET)
# ifndef __cray__
  /* Not sure why this function segfaults with Cray compiler */
  show_omp_thread_affinity(mpi);
# endif
#endif
}

void show_process_environment(int mpi)
{
#if defined(SPEC_OPENMP) || defined(SPEC_OPENMP_TARGET)
  char *omp_num_threads, *omp_places, *omp_proc_bind;
  omp_num_threads = getenv("OMP_NUM_THREADS");
  omp_places = getenv("OMP_PLACES");
  omp_proc_bind = getenv("OMP_PROC_BIND");
  fprintf(stderr, "mpi=%d OMP_NUM_THREADS=%s OMP_PLACES=%s OMP_PROC_BIND=%s\n",
	  mpi, omp_num_threads, omp_places, omp_proc_bind);
#endif
#if defined(SPEC_OPENMP_TARGET) || defined(SPEC_OPENACC) || defined(SPEC_CUDA)
  char *cuda_visible_devices;
  cuda_visible_devices = getenv("CUDA_VISIBLE_DEVICES");
  fprintf(stderr, "mpi=%d CUDA_VISIBLE_DEVICES=%s\n", mpi, cuda_visible_devices);
#endif
  show_process_affinity(mpi);
}

void show_process_affinity(int mpi)
{
  cpu_set_t coremask;
  char clbuf[7 * CPU_SETSIZE], hnbuf[64];

  memset(clbuf, 0, sizeof(clbuf));
  memset(hnbuf, 0, sizeof(hnbuf));
  (void)gethostname(hnbuf, sizeof(hnbuf));
  (void)sched_getaffinity(0, sizeof(coremask), &coremask);
  cpuset_to_cstr(&coremask, clbuf);
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
  int gpu, totgpu;
# ifdef SPEC_OPENACC
  acc_device_t my_device_type;
  my_device_type = acc_get_device_type();
  totgpu = acc_get_num_devices(my_device_type);
  gpu = acc_get_device_num(my_device_type);
# endif
# ifdef SPEC_OPENMP_TARGET
  totgpu = omp_get_num_devices();
  gpu = omp_get_default_device();
# endif
  fprintf(stderr, "mpi=%d cpu=%d gpu=%d totgpu=%d host=%s core_affinity=%s [MPI]\n",
	  mpi, sched_getcpu(), gpu, totgpu, hnbuf, clbuf);
#else
  fprintf(stderr, "mpi=%d cpu=%d host=%s core_affinity=%s [MPI]\n",
	  mpi, sched_getcpu(), hnbuf, clbuf);
#endif
}

#if defined(SPEC_OPENMP) || defined(SPEC_OPENMP_TARGET)
void show_omp_thread_affinity(int mpi)
{
  int tid;
  cpu_set_t coremask;
  char clbuf[7 * CPU_SETSIZE], hnbuf[64];

#pragma omp parallel private(tid, coremask, hnbuf, clbuf) shared(mpi)
  {
    tid = omp_get_thread_num();
    memset(clbuf, 0, sizeof(clbuf));
    memset(hnbuf, 0, sizeof(hnbuf));
    (void)gethostname(hnbuf, sizeof(hnbuf));
    (void)sched_getaffinity(0, sizeof(coremask), &coremask);
    cpuset_to_cstr(&coremask, clbuf);
# ifdef SPEC_OPENMP_TARGET
    fprintf(stderr, "mpi=%d tid=%d cpu=%d gpu=%d totgpu=%d host=%s core_affinity=%s [OpenMP]\n",
	    mpi, tid, sched_getcpu(), omp_get_default_device(), omp_get_num_devices(),
	    hnbuf, clbuf);
# else
    fprintf(stderr, "mpi=%d tid=%d cpu=%d host=%s core_affinity=%s [OpenMP]\n",
	    mpi, tid, sched_getcpu(), hnbuf, clbuf);
# endif
  }
}
#endif

/* Borrowed from util-linux-2.13-pre7/schedutils/taskset.c */
static char *cpuset_to_cstr(cpu_set_t *mask, char *str)
{
 char *ptr = str;
 int i, j, entry_made = 0;
 for (i = 0; i < CPU_SETSIZE; i++) {
   if (CPU_ISSET(i, mask)) {
     int run = 0;
     entry_made = 1;
     for (j = i + 1; j < CPU_SETSIZE; j++) {
       if (CPU_ISSET(j, mask)) run++;
       else break;
     }
     if (!run)
       sprintf(ptr, "%d,", i);
     else if (run == 1) {
       sprintf(ptr, "%d,%d,", i, i + 1);
       i++;
     } else {
       sprintf(ptr, "%d-%d,", i, i + run);
       i += run;
     }
     while (*ptr != 0) ptr++;
   }
 }
 ptr -= entry_made;
 *ptr = 0;
 return(str);
}

void show_procmem(void)
{
  double lmem, gmem;
  int rank;

  lmem = get_procmem();
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Allreduce(&lmem, &gmem, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  rank = 0;
  gmem = lmem;
#endif
  if (rank == 0) {
    fprintf(stderr, "\nTotal memory usage = %.3lf GiB\n", gmem);
  }
}

double get_procmem(void)
{
  FILE *fh;
  char proc_var[80];
  char *cp;
  long long unsigned int ibytes;
  double GiB = 0.0;

  fh = fopen("/proc/self/status","r");
  while(!feof(fh)) {
    fgets(proc_var,80,fh);
    cp = strstr(proc_var,"VmHWM:");
    if (cp) {sscanf(cp, "VmHWM:"" %llu",&ibytes );
      GiB = (double) ibytes / (1024.0 * 1024.0);
    }
  }
  fclose(fh);
  return GiB;
}
#endif
