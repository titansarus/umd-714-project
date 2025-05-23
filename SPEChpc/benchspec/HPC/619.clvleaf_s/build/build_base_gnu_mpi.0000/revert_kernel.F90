!Crown Copyright 2012 AWE.
!
! This file is part of CloverLeaf.
!
! CloverLeaf is free software: you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! CloverLeaf is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details.
!
! You should have received a copy of the GNU General Public License along with 
! CloverLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Fortran revert kernel.
!>  @author Wayne Gaudin
!>  @details Takes the half step field data used in the predictor and reverts
!>  it to the start of step data, ready for the corrector.
!>  Note that this does not seem necessary in this proxy-app but should be
!>  left in to remain relevant to the full method.

MODULE revert_kernel_module

CONTAINS

  SUBROUTINE revert_kernel(x_min,x_max,y_min,y_max,density0,density1,energy0,energy1)

    IMPLICIT NONE

    INTEGER :: x_min,x_max,y_min,y_max
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2)    :: density0
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2)    :: density1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2)    :: energy0
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2)    :: energy1

    INTEGER :: j,k

#ifdef SPEC_OPENMP_TARGET
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(tofrom: density0,energy0,density1,energy1)
#endif
#if defined(SPEC_OPENMP) 
    !$OMP PARALLEL
    !$OMP DO
#endif 
#ifdef SPEC_OPENACC
!$ACC DATA &
!$ACC PRESENT(density0,energy0,density1,energy1)
!$ACC KERNELS
!$ACC LOOP INDEPENDENT
#endif 
    DO k=y_min,y_max
#ifdef SPEC_OPENACC
!$ACC LOOP INDEPENDENT
#endif
      DO j=x_min,x_max
        density1(j,k)=density0(j,k)
        energy1(j,k)=energy0(j,k)
      ENDDO
    ENDDO
#ifdef SPEC_OPENACC
!$ACC END KERNELS
!$ACC END DATA
#endif

#if defined(SPEC_OPENMP) 
  !$OMP END DO
  !$OMP END PARALLEL
#endif 

  END SUBROUTINE revert_kernel

END MODULE revert_kernel_module
