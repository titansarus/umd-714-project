//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include "timers.h"

/* timers.h defines getTime as a macro function when using MPI or OpenMP */
#ifndef getTime
# if defined(SPEC_OPENACC) || (defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 200112L)
# include <stddef.h> /* defines NULL */
# include <sys/time.h> /* defines timeval */
double getTime() {
  struct timeval t;
  gettimeofday(&t,NULL);
  return t.tv_sec + (1.0E-6 * t.tv_usec);
}
# else
/* Stub implementation of getTime. We only need a valid implementation
   of getTime when compiling with -DVERBOSE */
double getTime() {
  return -1.0;
}
# endif
#endif
