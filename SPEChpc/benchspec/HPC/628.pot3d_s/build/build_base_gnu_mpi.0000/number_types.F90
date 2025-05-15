!#######################################################################
      module number_types
!
!-----------------------------------------------------------------------
! ****** Basic number types.
! ****** This module is used to set the default precision for REALs.
!-----------------------------------------------------------------------
!
      use iso_fortran_env
!
!-----------------------------------------------------------------------
!
      implicit none
!
      integer, parameter :: KIND_REAL_4=REAL32
      integer, parameter :: KIND_REAL_8=REAL64
      integer, parameter :: KIND_REAL_16=max(REAL128,REAL64)
!
      integer, parameter :: r_typ=KIND_REAL_8
!
      end module
!#######################################################################
