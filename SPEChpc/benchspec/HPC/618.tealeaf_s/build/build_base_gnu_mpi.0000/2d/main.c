#ifdef SPEC_OPENACC
  #define DIFFUSE_OVERLOAD
#endif

#ifdef SPEC_OPENMP_TARGET
  #define DIFFUSE_OVERLOAD
#endif

#include "comms.h"
#include "application.h"
#include "chunk.h"
#include "shared.h"
#include "drivers/drivers.h"
#ifdef SPEC
#include "specmpitime.h"
#endif
void settings_overload(Settings* settings, int argc, char** argv);

int main(int argc, char** argv)
{
  // Immediately initialise MPI
  initialise_comms(argc, argv);
SPEC_TIME_START
SPEC_TIME_START_INIT
  barrier();

  // Create the settings wrapper
  Settings* settings = (Settings*)malloc(sizeof(Settings));
  set_default_settings(settings);

  // Fill in rank information
  initialise_ranks(settings);

  barrier();

  // Perform initialisation steps
  Chunk* chunks;
  initialise_application(&chunks, settings);

  settings_overload(settings, argc, argv);
SPEC_TIME_STOP_INIT
SPEC_TIME_START_COMP
  // Perform the solve using default or overloaded diffuse
#ifndef DIFFUSE_OVERLOAD
  diffuse(chunks, settings);
#else
  diffuse_overload(chunks, settings);
#endif
SPEC_TIME_STOP_COMP

#ifndef SPEC
  // Print the kernel-level profiling results
  if(settings->rank == MASTER)
  {
    PRINT_PROFILING_RESULTS(settings->kernel_profile);
  }
#endif

  // Finalise the kernel
  kernel_finalise_driver(chunks, settings);

  // Finalise each individual chunk
  for(int cc = 0; cc < settings->num_chunks_per_rank; ++cc)
  {
    finalise_chunk(&(chunks[cc]));
    free(&(chunks[cc]));
  }

  // Finalise the application
  free(settings);
SPEC_TIME_STOP
SPEC_TIME_FINAL(true,0.0)
  finalise_comms();

  return EXIT_SUCCESS;
}

void settings_overload(Settings* settings, int argc, char** argv)
{
  for(int aa = 1; aa < argc; ++aa)
  {
    // Overload the solver
    if(strmatch(argv[aa], "-solver") ||
        strmatch(argv[aa], "--solver") || 
        strmatch(argv[aa], "-s"))
    {
      if(aa+1 == argc) break;
      if(strmatch(argv[aa+1], "cg")) settings->solver = CG_SOLVER;
      if(strmatch(argv[aa+1], "cheby")) settings->solver = CHEBY_SOLVER;
      if(strmatch(argv[aa+1], "ppcg")) settings->solver = PPCG_SOLVER;
      if(strmatch(argv[aa+1], "jacobi")) settings->solver = JACOBI_SOLVER;
    }
    else if(strmatch(argv[aa], "-x"))
    {
      if(aa+1 == argc) break;
      settings->grid_x_cells = atoi(argv[aa]);
    }
    else if(strmatch(argv[aa], "-y"))
    {
      if(aa+1 == argc) break;
      settings->grid_y_cells = atoi(argv[aa]);
    }
    else if(strmatch(argv[aa], "-help") ||
        strmatch(argv[aa], "--help") || 
        strmatch(argv[aa], "-h"))
    {
      print_and_log(settings, "tealeaf <options>\n");
      print_and_log(settings, "options:\n");
      print_and_log(settings, "\t-solver, --solver, -s:\n");
      print_and_log(settings, 
          "\t\tCan be 'cg', 'cheby', 'ppcg', or 'jacobi'\n");
      finalise_comms();
      exit(0);
    } 
  }
}
