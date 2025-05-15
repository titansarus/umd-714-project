#include <stdlib.h> /* defines getenv */
#include <sys/time.h> /* defines timeval */
#include <stdio.h> /* defines printf */

#define NINIT 64
int directives_runtime_init(void)
{
  double x[NINIT];
  double sum, expected_sum, diff, abs_diff;
  struct timeval start_time, stop_time;
  double elapsed_time;
  size_t i;
  int rtn, show_runtime;
#ifdef VERBOSE
  show_runtime = 1;
#else
  show_runtime = 0;
#endif

  for (i=0; i<NINIT; ++i) x[i] = (double)i;
  expected_sum = (NINIT/2.) * (NINIT - 1.);
  sum = 0.0;

  gettimeofday(&start_time,NULL);

#if defined(SPEC_OPENMP_TARGET)
# pragma omp target teams distribute parallel for map(to:x[0:NINIT]) map(tofrom:sum) reduction(+:sum)
#elif defined(SPEC_OPENACC)
# pragma acc parallel loop copyin(x[0:NINIT]) copy(sum) reduction(+:sum)
#endif
  for (i=0; i<NINIT; ++i) sum += x[i];

  gettimeofday(&stop_time,NULL);
  elapsed_time = (stop_time.tv_sec - start_time.tv_sec) + \
    (1.0E-6 * (stop_time.tv_usec - start_time.tv_usec));

  /* Calculate absolute value inline to avoid dependency
     on libm which causes compilation issues on x86
     when using Clang compiler and OpenMP target offload */
  diff = sum - expected_sum;
  abs_diff = diff >= 0.0 ? diff : -diff;
  rtn = abs_diff < 1.0E-6 ? 0 : 1;

  if (show_runtime) {
    fprintf(stderr, "Runtime initialization time: %.6lf seconds (rtn=%d)\n",
	    elapsed_time, rtn);
  }
  return rtn;
}
