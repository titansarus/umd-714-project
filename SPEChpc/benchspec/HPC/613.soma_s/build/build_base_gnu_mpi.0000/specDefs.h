#pragma once

#ifndef SPEC_NO_VAR_ARRAY_REDUCE
#if ( defined SPEC_OPENACC && ( _OPENACC < 201811 ) ) \
	|| ( ( defined SPEC_OPENMP_TARGET || defined SPEC_OPENMP ) \
		&& ( _OPENMP < 201511 ) )
// variable size array reduce not supported by compiler
#define SPEC_NO_VAR_ARRAY_REDUCE
#endif
#endif
