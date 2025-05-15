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

!>  @brief Main set up routine
!>  @author Wayne Gaudin
!>  @details Invokes the mesh decomposer and sets up chunk connectivity. It then
!>  allocates the communication buffers and call the chunk initialisation and
!>  generation routines. It calls the equation of state to calculate initial
!>  pressure before priming the halo cells and writing an initial field summary.

SUBROUTINE start

  USE clover_module
  USE parse_module
  USE update_halo_module
  USE ideal_gas_module

  IMPLICIT NONE

  INTEGER :: c, tile

  INTEGER :: x_cells,y_cells
  INTEGER:: right,left,top,bottom

  INTEGER :: fields(NUM_FIELDS) !, chunk_task_responsible_for

  LOGICAL :: profiler_off

  IF(parallel%boss)THEN
    WRITE(g_out,*) 'Setting up initial geometry'
    WRITE(g_out,*)
  ENDIF

  time  = 0.0
  step  = 0
  dtold = dtinit
  dt    = dtinit
  CALL clover_barrier

  CALL clover_get_num_chunks(number_of_chunks)

  CALL clover_decompose(grid%x_cells,grid%y_cells,left,right,bottom,top)

  !create the chunks
      
  chunk%task = parallel%task

  !chunk_task_responsible_for = parallel%task+1

  x_cells = right -left  +1
  y_cells = top   -bottom+1
      
  chunk%left    = left
  chunk%bottom  = bottom
  chunk%right   = right
  chunk%top     = top
  chunk%left_boundary   = 1
  chunk%bottom_boundary = 1
  chunk%right_boundary  = grid%x_cells
  chunk%top_boundary    = grid%y_cells
  chunk%x_min = 1
  chunk%y_min = 1
  chunk%x_max = x_cells
  chunk%y_max = y_cells
    
  ! create the tiles
  ALLOCATE( chunk%tiles(1:tiles_per_chunk) )

  CALL clover_tile_decompose(x_cells, y_cells)

  CALL build_field()

  CALL clover_barrier
  CALL clover_allocate_buffers()

  IF(parallel%boss)THEN
    WRITE(g_out,*) 'Generating chunks'
  ENDIF

#ifdef SPEC_OPENACC
!$ACC enter DATA &
    !$ACC CREATE(chunk%tiles(1)%field%density0)   &
    !$ACC CREATE(chunk%tiles(1)%field%density1)   &
    !$ACC CREATE(chunk%tiles(1)%field%energy0)    &
    !$ACC CREATE(chunk%tiles(1)%field%energy1)    &
    !$ACC CREATE(chunk%tiles(1)%field%pressure)   &
    !$ACC CREATE(chunk%tiles(1)%field%soundspeed) &
    !$ACC CREATE(chunk%tiles(1)%field%viscosity)  &
    !$ACC CREATE(chunk%tiles(1)%field%xvel0)      &
    !$ACC CREATE(chunk%tiles(1)%field%yvel0)      &
    !$ACC CREATE(chunk%tiles(1)%field%xvel1)      &
    !$ACC CREATE(chunk%tiles(1)%field%yvel1)      &
    !$ACC CREATE(chunk%tiles(1)%field%vol_flux_x) &
    !$ACC CREATE(chunk%tiles(1)%field%vol_flux_y) &
    !$ACC CREATE(chunk%tiles(1)%field%mass_flux_x)&
    !$ACC CREATE(chunk%tiles(1)%field%mass_flux_y)&
    !$ACC CREATE(chunk%tiles(1)%field%volume)     &
    !$ACC CREATE(chunk%tiles(1)%field%work_array1)&
    !$ACC CREATE(chunk%tiles(1)%field%work_array2)&
    !$ACC CREATE(chunk%tiles(1)%field%work_array3)&
    !$ACC CREATE(chunk%tiles(1)%field%work_array4)&
    !$ACC CREATE(chunk%tiles(1)%field%work_array5)&
    !$ACC CREATE(chunk%tiles(1)%field%work_array6)&
    !$ACC CREATE(chunk%tiles(1)%field%work_array7)&
    !$ACC CREATE(chunk%tiles(1)%field%cellx)      &
    !$ACC CREATE(chunk%tiles(1)%field%celly)      &
    !$ACC CREATE(chunk%tiles(1)%field%celldx)     &
    !$ACC CREATE(chunk%tiles(1)%field%celldy)     &
    !$ACC CREATE(chunk%tiles(1)%field%vertexx)    &
    !$ACC CREATE(chunk%tiles(1)%field%vertexdx)   &
    !$ACC CREATE(chunk%tiles(1)%field%vertexy)    &
    !$ACC CREATE(chunk%tiles(1)%field%vertexdy)   &
    !$ACC CREATE(chunk%tiles(1)%field%xarea)      &
    !$ACC CREATE(chunk%tiles(1)%field%yarea)      &
    !$ACC CREATE(chunk%left_snd_buffer)    &
    !$ACC CREATE(chunk%left_rcv_buffer)    &
    !$ACC CREATE(chunk%right_snd_buffer)   &
    !$ACC CREATE(chunk%right_rcv_buffer)   &
    !$ACC CREATE(chunk%bottom_snd_buffer)  &
    !$ACC CREATE(chunk%bottom_rcv_buffer)  &
    !$ACC CREATE(chunk%top_snd_buffer)     &
    !$ACC CREATE(chunk%top_rcv_buffer)
#endif
#ifdef SPEC_OPENMP_TARGET
!$omp target enter data map(alloc:       &
!$omp   chunk%tiles(1)%field%density0,   &
!$omp   chunk%tiles(1)%field%density1,   &
!$omp   chunk%tiles(1)%field%energy0,    &
!$omp   chunk%tiles(1)%field%energy1,    &
!$omp   chunk%tiles(1)%field%pressure,   &
!$omp   chunk%tiles(1)%field%soundspeed, &
!$omp   chunk%tiles(1)%field%viscosity,  &
!$omp   chunk%tiles(1)%field%xvel0,      &
!$omp   chunk%tiles(1)%field%yvel0,      &
!$omp   chunk%tiles(1)%field%xvel1,      &
!$omp   chunk%tiles(1)%field%yvel1,      &
!$omp   chunk%tiles(1)%field%vol_flux_x, &
!$omp   chunk%tiles(1)%field%vol_flux_y, &
!$omp   chunk%tiles(1)%field%mass_flux_x,&
!$omp   chunk%tiles(1)%field%mass_flux_y,&
!$omp   chunk%tiles(1)%field%volume,     &
!$omp   chunk%tiles(1)%field%work_array1,&
!$omp   chunk%tiles(1)%field%work_array2,&
!$omp   chunk%tiles(1)%field%work_array3,&
!$omp   chunk%tiles(1)%field%work_array4,&
!$omp   chunk%tiles(1)%field%work_array5,&
!$omp   chunk%tiles(1)%field%work_array6,&
!$omp   chunk%tiles(1)%field%work_array7,&
!$omp   chunk%tiles(1)%field%cellx,      &
!$omp   chunk%tiles(1)%field%celly,      &
!$omp   chunk%tiles(1)%field%celldx,     &
!$omp   chunk%tiles(1)%field%celldy,     &
!$omp   chunk%tiles(1)%field%vertexx,    &
!$omp   chunk%tiles(1)%field%vertexdx,   &
!$omp   chunk%tiles(1)%field%vertexy,    &
!$omp   chunk%tiles(1)%field%vertexdy,   &
!$omp   chunk%tiles(1)%field%xarea,      &
!$omp   chunk%tiles(1)%field%yarea,      &
!$omp   chunk%left_snd_buffer,    &
!$omp   chunk%left_rcv_buffer,    &
!$omp   chunk%right_snd_buffer,   &
!$omp   chunk%right_rcv_buffer,   &
!$omp   chunk%bottom_snd_buffer,  &
!$omp   chunk%bottom_rcv_buffer,  &
!$omp   chunk%top_snd_buffer,     &
!$omp   chunk%top_rcv_buffer)
#endif


  DO tile=1,tiles_per_chunk
    CALL initialise_chunk(tile)
    CALL generate_chunk(tile)
  ENDDO

  advect_x=.TRUE.
  CALL clover_barrier

  ! Do no profile the start up costs otherwise the total times will not add up
  ! at the end
  profiler_off=profiler_on
  profiler_on=.FALSE.

  DO tile = 1, tiles_per_chunk
    CALL ideal_gas(tile,.FALSE.)
  END DO

  ! Prime all halo data for the first step
  fields=0
  fields(FIELD_DENSITY0)=1
  fields(FIELD_ENERGY0)=1
  fields(FIELD_PRESSURE)=1
  fields(FIELD_VISCOSITY)=1
  fields(FIELD_DENSITY1)=1
  fields(FIELD_ENERGY1)=1
  fields(FIELD_XVEL0)=1
  fields(FIELD_YVEL0)=1
  fields(FIELD_XVEL1)=1
  fields(FIELD_YVEL1)=1
  CALL update_halo(fields,2)

  IF(parallel%boss)THEN
    WRITE(g_out,*)
    WRITE(g_out,*) 'Problem initialised and generated'
  ENDIF
  CALL field_summary()

  IF(visit_frequency.NE.0) CALL visit()

#ifdef SPEC_OPENACC
!$ACC exit DATA &
    !$ACC COPYOUT(chunk%tiles(1)%field%density0)   &
    !$ACC COPYOUT(chunk%tiles(1)%field%density1)   &
    !$ACC COPYOUT(chunk%tiles(1)%field%energy0)    &
    !$ACC COPYOUT(chunk%tiles(1)%field%energy1)    &
    !$ACC COPYOUT(chunk%tiles(1)%field%pressure)   &
    !$ACC COPYOUT(chunk%tiles(1)%field%soundspeed) &
    !$ACC COPYOUT(chunk%tiles(1)%field%viscosity)  &
    !$ACC COPYOUT(chunk%tiles(1)%field%xvel0)      &
    !$ACC COPYOUT(chunk%tiles(1)%field%yvel0)      &
    !$ACC COPYOUT(chunk%tiles(1)%field%xvel1)      &
    !$ACC COPYOUT(chunk%tiles(1)%field%yvel1)      &
    !$ACC COPYOUT(chunk%tiles(1)%field%vol_flux_x) &
    !$ACC COPYOUT(chunk%tiles(1)%field%vol_flux_y) &
    !$ACC COPYOUT(chunk%tiles(1)%field%mass_flux_x)&
    !$ACC COPYOUT(chunk%tiles(1)%field%mass_flux_y)&
    !$ACC COPYOUT(chunk%tiles(1)%field%volume)     &
    !$ACC COPYOUT(chunk%tiles(1)%field%work_array1)&
    !$ACC COPYOUT(chunk%tiles(1)%field%work_array2)&
    !$ACC COPYOUT(chunk%tiles(1)%field%work_array3)&
    !$ACC COPYOUT(chunk%tiles(1)%field%work_array4)&
    !$ACC COPYOUT(chunk%tiles(1)%field%work_array5)&
    !$ACC COPYOUT(chunk%tiles(1)%field%work_array6)&
    !$ACC COPYOUT(chunk%tiles(1)%field%work_array7)&
    !$ACC COPYOUT(chunk%tiles(1)%field%cellx)      &
    !$ACC COPYOUT(chunk%tiles(1)%field%celly)      &
    !$ACC COPYOUT(chunk%tiles(1)%field%celldx)     &
    !$ACC COPYOUT(chunk%tiles(1)%field%celldy)     &
    !$ACC COPYOUT(chunk%tiles(1)%field%vertexx)    &
    !$ACC COPYOUT(chunk%tiles(1)%field%vertexdx)   &
    !$ACC COPYOUT(chunk%tiles(1)%field%vertexy)    &
    !$ACC COPYOUT(chunk%tiles(1)%field%vertexdy)   &
    !$ACC COPYOUT(chunk%tiles(1)%field%xarea)      &
    !$ACC COPYOUT(chunk%tiles(1)%field%yarea)      &
    !$ACC COPYOUT(chunk%left_snd_buffer)    &
    !$ACC COPYOUT(chunk%left_rcv_buffer)    &
    !$ACC COPYOUT(chunk%right_snd_buffer)   &
    !$ACC COPYOUT(chunk%right_rcv_buffer)   &
    !$ACC COPYOUT(chunk%bottom_snd_buffer)  &
    !$ACC COPYOUT(chunk%bottom_rcv_buffer)  &
    !$ACC COPYOUT(chunk%top_snd_buffer)     &
    !$ACC COPYOUT(chunk%top_rcv_buffer)
#endif   
#ifdef SPEC_OPENMP_TARGET
!$omp target exit data map(from:         &
!$omp   chunk%tiles(1)%field%density0,   &
!$omp   chunk%tiles(1)%field%density1,   &
!$omp   chunk%tiles(1)%field%energy0,    &
!$omp   chunk%tiles(1)%field%energy1,    &
!$omp   chunk%tiles(1)%field%pressure,   &
!$omp   chunk%tiles(1)%field%soundspeed, &
!$omp   chunk%tiles(1)%field%viscosity,  &
!$omp   chunk%tiles(1)%field%xvel0,      &
!$omp   chunk%tiles(1)%field%yvel0,      &
!$omp   chunk%tiles(1)%field%xvel1,      &
!$omp   chunk%tiles(1)%field%yvel1,      &
!$omp   chunk%tiles(1)%field%vol_flux_x, &
!$omp   chunk%tiles(1)%field%vol_flux_y, &
!$omp   chunk%tiles(1)%field%mass_flux_x,&
!$omp   chunk%tiles(1)%field%mass_flux_y,&
!$omp   chunk%tiles(1)%field%volume,     &
!$omp   chunk%tiles(1)%field%work_array1,&
!$omp   chunk%tiles(1)%field%work_array2,&
!$omp   chunk%tiles(1)%field%work_array3,&
!$omp   chunk%tiles(1)%field%work_array4,&
!$omp   chunk%tiles(1)%field%work_array5,&
!$omp   chunk%tiles(1)%field%work_array6,&
!$omp   chunk%tiles(1)%field%work_array7,&
!$omp   chunk%tiles(1)%field%cellx,      &
!$omp   chunk%tiles(1)%field%celly,      &
!$omp   chunk%tiles(1)%field%celldx,     &
!$omp   chunk%tiles(1)%field%celldy,     &
!$omp   chunk%tiles(1)%field%vertexx,    &
!$omp   chunk%tiles(1)%field%vertexdx,   &
!$omp   chunk%tiles(1)%field%vertexy,    &
!$omp   chunk%tiles(1)%field%vertexdy,   &
!$omp   chunk%tiles(1)%field%xarea,      &
!$omp   chunk%tiles(1)%field%yarea,      &
!$omp   chunk%left_snd_buffer,    &
!$omp   chunk%left_rcv_buffer,    &
!$omp   chunk%right_snd_buffer,   &
!$omp   chunk%right_rcv_buffer,   &
!$omp   chunk%bottom_snd_buffer,  &
!$omp   chunk%bottom_rcv_buffer,  &
!$omp   chunk%top_snd_buffer,     &
!$omp   chunk%top_rcv_buffer)

#endif
  CALL clover_barrier

  profiler_on=profiler_off

END SUBROUTINE start
