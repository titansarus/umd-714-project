/* Copyright (C) 2016 Ludwig Schneider

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

//! \file soma_config.h
//! \brief configuration variables for SOMA File is configured by CMAKE.
#ifndef SOMA_CONFIG_H
#define SOMA_CONFIG_H

//! Macro to destinguish between SINGLE_PRECISION and DOUBLE.
//! Automatically set by CMake
//#cmakedefine01 SINGLE_PRECISION
#define SINGLE_PRECISION 0

#if ( SINGLE_PRECISION == 1)

//! Alias for float or double variables
typedef float soma_scalar_t;
//! Alias for float or double variables in HDF5 memory
#define H5T_SOMA_NATIVE_SCALAR H5T_NATIVE_FLOAT
//! Alias for float or double variables in HDF5 files
#define H5T_SOMA_FILE_SCALAR H5T_IEEE_F32LE
//! Alias for flow or double variable in MPI
#define MPI_SOMA_SCALAR MPI_FLOAT

#else//SINGLE_PRECISION

//! Alias for float or double variables
typedef double soma_scalar_t;
//! Alias for float of double variables in HDF5 memory
#define H5T_SOMA_NATIVE_SCALAR H5T_NATIVE_DOUBLE
//! Alias for float of double variables in HDF5 files
#define H5T_SOMA_FILE_SCALAR H5T_IEEE_F64LE
//! Alias for flow or double variable in MPI
#define MPI_SOMA_SCALAR MPI_DOUBLE

#endif//SINGLE_PRECISION

//! String containing the git version of SOMA
static const char soma_version[] = "benchmark version based on tag 0.5.0 gitlab.com/InnocentBug/SOMA";
//! Returns the version string.
//! \return Pointer to version string
const char *get_soma_version(void);

//! String describing the system info for which SOMA has been compiled.
static const char soma_system_info[] = "system info not available";
//! Returns the string describing the system SOMA has been compiled.
//! \return Pointer to string.
const char *get_soma_system_info(void);

#endif//SOMA_CONFIG_H

//Check for OMP or OPENACC usage
#ifdef SPEC_OPENACC
#ifdef SPEC_OPENMP
#error "You could either compile with OPENACC or OMP. Not with both."
#endif//SPEC_OPENMP
#endif//SPEC_OPENACC
