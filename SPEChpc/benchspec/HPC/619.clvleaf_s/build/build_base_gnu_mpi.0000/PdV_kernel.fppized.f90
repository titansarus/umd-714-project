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

!>  @brief Fortran PdV kernel.
!>  @author Wayne Gaudin
!>  @details Calculates the change in energy and density in a cell using the
!>  change on cell volume due to the velocity gradients in a cell. The time
!>  level of the velocity data depends on whether it is invoked as the
!>  predictor or corrector.

MODULE PdV_kernel_module

CONTAINS

  SUBROUTINE PdV_kernel(predict,                                          &
                        x_min,x_max,y_min,y_max,dt,                       &
                        xarea,yarea,volume,                               &
                        density0,                                         &
                        density1,                                         &
                        energy0,                                          &
                        energy1,                                          &
                        pressure,                                         &
                        viscosity,                                        &
                        xvel0,                                            &
                        xvel1,                                            &
                        yvel0,                                            &
                        yvel1,                                            &
                        volume_change                                     )

    IMPLICIT NONE

    LOGICAL :: predict

    INTEGER :: x_min,x_max,y_min,y_max
    REAL(KIND=8)  :: dt
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: xarea
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: yarea
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: volume
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0,energy0
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: pressure
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density1,energy1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: viscosity
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0,yvel0
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel1,yvel1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: volume_change

    INTEGER :: j,k

    REAL(KIND=8)  :: recip_volume,energy_change,min_cell_volume
    REAL(KIND=8)  :: right_flux,left_flux,top_flux,bottom_flux,total_flux
    REAL(KIND=8)  :: volume_change_s



    IF(predict)THEN

      DO k=y_min,y_max
        DO j=x_min,x_max

          left_flux=  (xarea(j  ,k  )*(xvel0(j  ,k  )+xvel0(j  ,k+1)                     &
            +xvel0(j  ,k  )+xvel0(j  ,k+1)))*0.25_8*dt*0.5
          right_flux= (xarea(j+1,k  )*(xvel0(j+1,k  )+xvel0(j+1,k+1)                     &
            +xvel0(j+1,k  )+xvel0(j+1,k+1)))*0.25_8*dt*0.5
          bottom_flux=(yarea(j  ,k  )*(yvel0(j  ,k  )+yvel0(j+1,k  )                     &
            +yvel0(j  ,k  )+yvel0(j+1,k  )))*0.25_8*dt*0.5
          top_flux=   (yarea(j  ,k+1)*(yvel0(j  ,k+1)+yvel0(j+1,k+1)                     &
            +yvel0(j  ,k+1)+yvel0(j+1,k+1)))*0.25_8*dt*0.5
          total_flux=right_flux-left_flux+top_flux-bottom_flux

          volume_change_s=volume(j,k)/(volume(j,k)+total_flux)

          min_cell_volume=MIN(volume(j,k)+right_flux-left_flux+top_flux-bottom_flux &
            ,volume(j,k)+right_flux-left_flux                      &
            ,volume(j,k)+top_flux-bottom_flux)
 
          recip_volume=1.0/volume(j,k)

          energy_change=(pressure(j,k)/density0(j,k)+viscosity(j,k)/density0(j,k))*total_flux*recip_volume

          energy1(j,k)=energy0(j,k)-energy_change

          density1(j,k)=density0(j,k)*volume_change_s

        ENDDO
      ENDDO



    ELSE


      DO k=y_min,y_max
        DO j=x_min,x_max

          left_flux=  (xarea(j  ,k  )*(xvel0(j  ,k  )+xvel0(j  ,k+1)                     &
            +xvel1(j  ,k  )+xvel1(j  ,k+1)))*0.25_8*dt
          right_flux= (xarea(j+1,k  )*(xvel0(j+1,k  )+xvel0(j+1,k+1)                     &
            +xvel1(j+1,k  )+xvel1(j+1,k+1)))*0.25_8*dt
          bottom_flux=(yarea(j  ,k  )*(yvel0(j  ,k  )+yvel0(j+1,k  )                     &
            +yvel1(j  ,k  )+yvel1(j+1,k  )))*0.25_8*dt
          top_flux=   (yarea(j  ,k+1)*(yvel0(j  ,k+1)+yvel0(j+1,k+1)                     &
            +yvel1(j  ,k+1)+yvel1(j+1,k+1)))*0.25_8*dt
          total_flux=right_flux-left_flux+top_flux-bottom_flux

          volume_change_s=volume(j,k)/(volume(j,k)+total_flux)

          min_cell_volume=MIN(volume(j,k)+right_flux-left_flux+top_flux-bottom_flux &
            ,volume(j,k)+right_flux-left_flux                      &
            ,volume(j,k)+top_flux-bottom_flux)
 
          recip_volume=1.0/volume(j,k)

          energy_change=(pressure(j,k)/density0(j,k)+viscosity(j,k)/density0(j,k))*total_flux*recip_volume

          energy1(j,k)=energy0(j,k)-energy_change

          density1(j,k)=density0(j,k)*volume_change_s

        ENDDO
      ENDDO


    ENDIF




  END SUBROUTINE PdV_kernel

END MODULE PdV_kernel_module

