! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!   Copyright by The HDF Group.                                               *
!   All rights reserved.                                                      *
!                                                                             *
!   This file is part of HDF5.  The full HDF5 copyright notice, including     *
!   terms governing use, modification, and redistribution, is contained in    *
!   the COPYING file, which can be found at the root of the source code       *
!   distribution tree, or in https://support.hdfgroup.org/ftp/HDF5/releases.  *
!   If you do not have access to either file, you may request a copy from     *
!   help@hdfgroup.org.                                                        *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! fortran/src/H5config_f.inc. Generated from fortran/src/H5config_f.inc.in by configure

! The script to replace the defines in H5config_f.inc.in is
! located in configure.ac in the Fortran section. 

! Define if we have parallel support
! SPEC default is enabled
#ifndef SPEC_NO_H5_HAVE_PARALLEL
#define H5_HAVE_PARALLEL 1
#else
#undef H5_HAVE_PARALLEL
#endif

! Define if the intrinsic function STORAGE_SIZE exists
! SPEC default is enabled
#ifndef SPEC_NO_H5_FORTRAN_HAVE_STORAGE_SIZE
#define H5_FORTRAN_HAVE_STORAGE_SIZE 1
#else
#undef H5_FORTRAN_HAVE_STORAGE_SIZE
#endif

! Define if the intrinsic function SIZEOF exists
! SPEC default is enabled
#ifndef SPEC_NO_H5_FORTRAN_HAVE_SIZEOF
#define H5_FORTRAN_HAVE_SIZEOF 1
#else
#undef H5_FORTRAN_HAVE_SIZEOF
#endif

! Define if the intrinsic function C_SIZEOF exists
! SPEC default is enabled
#ifndef SPEC_NO_H5_FORTRAN_HAVE_C_SIZEOF
#define H5_FORTRAN_HAVE_C_SIZEOF 1
#else
#undef H5_FORTRAN_HAVE_C_SIZEOF
#endif

! Define if the intrinsic function C_LONG_DOUBLE exists
! SPEC default is enabled
#ifndef SPEC_NO_H5_FORTRAN_HAVE_C_LONG_DOUBLE
#define H5_FORTRAN_HAVE_C_LONG_DOUBLE 1
#else
#undef H5_FORTRAN_HAVE_C_LONG_DOUBLE
#endif

! Define if Fortran C_LONG_DOUBLE is different from C_DOUBLE
! SPEC default is disabled 
#ifdef SPEC_H5_FORTRAN_C_LONG_DOUBLE_IS_UNIQUE
#define H5_FORTRAN_C_LONG_DOUBLE_IS_UNIQUE SPEC_H5_FORTRAN_C_LONG_DOUBLE_IS_UNIQUE
#else
#define H5_FORTRAN_C_LONG_DOUBLE_IS_UNIQUE 0
#endif

! Define if the intrinsic module ISO_FORTRAN_ENV exists
! SPEC default is disabled 
#ifdef SPEC_H5_HAVE_ISO_FORTRAN_ENV
#define H5_HAVE_ISO_FORTRAN_ENV 1
#else
#undef H5_HAVE_ISO_FORTRAN_ENV
#endif

! Define the size of C's double
! SPEC default is for doubles == 8
#ifdef SPEC_H5_SIZEOF_DOUBLE
#define H5_SIZEOF_DOUBLE SPEC_H5_SIZEOF_DOUBLE
#else
#define H5_SIZEOF_DOUBLE 8
#endif

! Define the size of C's long double
! SPEC default is for long doubles == 16
#ifdef SPEC_H5_SIZEOF_LONG_DOUBLE
#define H5_SIZEOF_LONG_DOUBLE SPEC_H5_SIZEOF_LONG_DOUBLE
#else
#define H5_SIZEOF_LONG_DOUBLE 16
#endif

! Define the maximum decimal precision for reals
#define H5_PAC_FC_MAX_REAL_PRECISION 15

! If C has quad precision
! SPEC default is disabled
#ifdef SPEC_H5_HAVE_FLOAT128
#define H5_HAVE_FLOAT128 1
#else
#undef H5_HAVE_FLOAT128
#endif

! Define if INTEGER*16 is available 
! SPEC default is disabled
#ifdef SPEC_H5_HAVE_Fortran_INTEGER_SIZEOF_16
#define H5_HAVE_Fortran_INTEGER_SIZEOF_16 1
#else
#define H5_HAVE_Fortran_INTEGER_SIZEOF_16 0
#endif

! Maximum decimal precision for C
! SPEC default is 21
#ifdef SPEC_H5_PAC_C_MAX_REAL_PRECISION
#define H5_PAC_C_MAX_REAL_PRECISION SPEC_H5_PAC_C_MAX_REAL_PRECISION
#else
#define H5_PAC_C_MAX_REAL_PRECISION 21
#endif

! number of valid REAL KINDs
#define H5_H5CONFIG_F_NUM_RKIND INTEGER, PARAMETER :: num_rkinds = 2

! valid REAL KINDs (need to have a matching C counter-part)
#define H5_H5CONFIG_F_RKIND INTEGER, DIMENSION(1:num_rkinds) :: rkind = (/4,8/)

! valid REAL KINDs (need to have a matching C counter-part)
#define H5_H5CONFIG_F_RKIND_SIZEOF INTEGER, DIMENSION(1:num_rkinds) :: rkind_sizeof = (/4,8/)

! number of valid INTEGER KINDs
#define H5_H5CONFIG_F_NUM_IKIND INTEGER, PARAMETER :: num_ikinds = 4

! valid INTEGER KINDs (need to have a matching C counter-part)
#define H5_H5CONFIG_F_IKIND INTEGER, DIMENSION(1:num_ikinds) :: ikind = (/1,2,4,8/)

! Fortran compiler name
#define H5_Fortran_COMPILER_ID none

