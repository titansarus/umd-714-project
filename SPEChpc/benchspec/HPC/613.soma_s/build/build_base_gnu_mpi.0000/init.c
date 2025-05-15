/* Copyright (C) 2016-2017 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016-2017 Marcel Langenberg
   Copyright (C) 2016 Fabien Leonforte
   Copyright (C) 2016 Juan Orozco
   Copyright (C) 2016 Yongzhi Ren
   Copyright (C) 2016 N. Harshavardhan Reddy

 This file is part of SOMA.

 SOMA is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SOMA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SOMA.  If not, see <http://www.gnu.org/licenses/>.
*/

//! \file init.c
//! \brief Implementation of init.h

#include "init.h"
#include <math.h>
#include <time.h>
#include <assert.h>
#include <memory.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <mpi.h>
#include <stdio.h>
#ifdef SPEC_OPENACC
#include <openacc.h>
#endif//SPEC_OPENACC
#if (defined(SPEC_OPENMP) || defined(SPEC_OPENMP_TARGET))
#include <omp.h>
#endif//SPEC_OPENMP
#include "soma_config.h"
#include "phase.h"

int print_version(const int rank)
    {
    if (rank == 0)
	{
	//Sytem
	fprintf(stdout,"system is %s with C std %ld.\n",get_soma_system_info(),__STDC_VERSION__);

	//SOMA git
	fprintf(stdout, "GIT version of SOMA is %s compiled on %s %s.\n",get_soma_version(),__DATE__,__TIME__);

#ifdef MPI_MAX_LIBRARY_VERSION_STRING
	//MPI
	char mpi_version[MPI_MAX_LIBRARY_VERSION_STRING];
	int length;
	MPI_Get_library_version(mpi_version,&length);
	int mpi_maj=0,mpi_min=0;
	MPI_Get_version(&mpi_maj,&mpi_min);
	fprintf(stdout,"MPI version: %s %d.%d\n",mpi_version,mpi_maj,mpi_min);
#else
	fprintf(stdout,"No MPI lib version available.\n");
#endif//mpi_max_library_version_string
	}
    return 0;
    }
