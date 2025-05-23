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

!>  @brief Controls the main hydro cycle.
!>  @author Wayne Gaudin
!>  @details Controls the top level cycle, invoking all the drivers and checks
!>  for outputs and completion.

SUBROUTINE hydro

    USE clover_module
    USE timestep_module
    USE viscosity_module
    USE PdV_module
    USE accelerate_module
    USE flux_calc_module
    USE advection_module
    USE reset_field_module

    IMPLICIT NONE

    INTEGER         :: loc(1)
    REAL(KIND=8)    :: timer,timerstart,wall_clock,step_clock
  
    REAL(KIND=8)    :: grind_time,cells,rstep
    REAL(KIND=8)    :: step_time,step_grind
    REAL(KIND=8)    :: first_step,second_step
    REAL(KIND=8)    :: kernel_total,totals(parallel%max_task)

    timerstart = timer()
    DO

        step_time = timer()

        step = step + 1
        CALL timestep()

        CALL PdV(.TRUE.)

        CALL accelerate()

        CALL PdV(.FALSE.)

        CALL flux_calc()

        CALL advection()

        CALL reset_field()

        advect_x = .NOT. advect_x
  
        time = time + dt

        IF(summary_frequency.NE.0) THEN
            IF(MOD(step, summary_frequency).EQ.0) CALL field_summary()
        ENDIF
        IF(visit_frequency.NE.0) THEN
            IF(MOD(step, visit_frequency).EQ.0) CALL visit()
        ENDIF

        ! Sometimes there can be a significant start up cost that appears in the first step.
        ! Sometimes it is due to the number of MPI tasks, or OpenCL kernel compilation.
        ! On the short test runs, this can skew the results, so should be taken into account
        !  in recorded run times.
        IF(step.EQ.1) first_step=(timer() - step_time)
        IF(step.EQ.2) second_step=(timer() - step_time)

        IF(time+g_small.GT.end_time.OR.step.GE.end_step) THEN

            complete=.TRUE.
            CALL field_summary()
            IF(visit_frequency.NE.0) CALL visit()


            IF ( profiler_on ) THEN

                ! First we need to find the maximum kernel time for each task. This
                ! seems to work better than finding the maximum time for each kernel and
                ! adding it up, which always gives over 100%. I think this is because it
                ! does not take into account compute overlaps before syncronisations
                ! caused by halo exhanges.
                kernel_total=profiler%timestep+profiler%ideal_gas+profiler%viscosity+profiler%PdV          &
                    +profiler%revert+profiler%acceleration+profiler%flux+profiler%cell_advection   &
                    +profiler%mom_advection+profiler%reset+profiler%summary+profiler%visit         &
                    +profiler%tile_halo_exchange+profiler%self_halo_exchange+profiler%mpi_halo_exchange
                CALL clover_allgather(kernel_total,totals)
                ! So then what I do is use the individual kernel times for the
                ! maximum kernel time task for the profile print
                loc=MAXLOC(totals)
                kernel_total=totals(loc(1))
                CALL clover_allgather(profiler%timestep,totals)
                profiler%timestep=totals(loc(1))
                CALL clover_allgather(profiler%ideal_gas,totals)
                profiler%ideal_gas=totals(loc(1))
                CALL clover_allgather(profiler%viscosity,totals)
                profiler%viscosity=totals(loc(1))
                CALL clover_allgather(profiler%PdV,totals)
                profiler%PdV=totals(loc(1))
                CALL clover_allgather(profiler%revert,totals)
                profiler%revert=totals(loc(1))
                CALL clover_allgather(profiler%acceleration,totals)
                profiler%acceleration=totals(loc(1))
                CALL clover_allgather(profiler%flux,totals)
                profiler%flux=totals(loc(1))
                CALL clover_allgather(profiler%cell_advection,totals)
                profiler%cell_advection=totals(loc(1))
                CALL clover_allgather(profiler%mom_advection,totals)
                profiler%mom_advection=totals(loc(1))
                CALL clover_allgather(profiler%reset,totals)
                profiler%reset=totals(loc(1))
                CALL clover_allgather(profiler%tile_halo_exchange,totals)
                profiler%tile_halo_exchange=totals(loc(1))
                CALL clover_allgather(profiler%self_halo_exchange,totals)
                profiler%self_halo_exchange=totals(loc(1))
                CALL clover_allgather(profiler%mpi_halo_exchange,totals)
                profiler%mpi_halo_exchange=totals(loc(1))
                CALL clover_allgather(profiler%summary,totals)
                profiler%summary=totals(loc(1))
                CALL clover_allgather(profiler%visit,totals)
                profiler%visit=totals(loc(1))

                IF ( parallel%boss ) THEN
                    WRITE(g_out,*)
                    WRITE(g_out,'(a58,2f16.4)')"Profiler Output                 Time            Percentage"
                    WRITE(g_out,'(a23,2f16.4)')"Timestep              :",profiler%timestep,&
                        100.0*(profiler%timestep/wall_clock)
                    WRITE(g_out,'(a23,2f16.4)')"Ideal Gas             :",profiler%ideal_gas,&
                        100.0*(profiler%ideal_gas/wall_clock)
                    WRITE(g_out,'(a23,2f16.4)')"Viscosity             :",profiler%viscosity,&
                        100.0*(profiler%viscosity/wall_clock)
                    WRITE(g_out,'(a23,2f16.4)')"PdV                   :",profiler%PdV,&
                        100.0*(profiler%PdV/wall_clock)
                    WRITE(g_out,'(a23,2f16.4)')"Revert                :",profiler%revert,&
                        100.0*(profiler%revert/wall_clock)
                    WRITE(g_out,'(a23,2f16.4)')"Acceleration          :",profiler%acceleration,&
                        100.0*(profiler%acceleration/wall_clock)
                    WRITE(g_out,'(a23,2f16.4)')"Fluxes                :",profiler%flux,&
                        100.0*(profiler%flux/wall_clock)
                    WRITE(g_out,'(a23,2f16.4)')"Cell Advection        :",profiler%cell_advection,&
                        100.0*(profiler%cell_advection/wall_clock)
                    WRITE(g_out,'(a23,2f16.4)')"Momentum Advection    :",profiler%mom_advection,&
                        100.0*(profiler%mom_advection/wall_clock)
                    WRITE(g_out,'(a23,2f16.4)')"Reset                 :",profiler%reset,&
                        100.0*(profiler%reset/wall_clock)
                    WRITE(g_out,'(a23,2f16.4)')"Summary               :",profiler%summary,&
                        100.0*(profiler%summary/wall_clock)
                    WRITE(g_out,'(a23,2f16.4)')"Visit                 :",profiler%visit,&
                        100.0*(profiler%visit/wall_clock)
                    WRITE(g_out,'(a23,2f16.4)')"Tile Halo Exchange    :",profiler%tile_halo_exchange,&
                        100.0*(profiler%tile_halo_exchange/wall_clock)
                    WRITE(g_out,'(a23,2f16.4)')"Self Halo Exchange    :",profiler%self_halo_exchange,&
                        100.0*(profiler%self_halo_exchange/wall_clock)
                    WRITE(g_out,'(a23,2f16.4)')"MPI Halo Exchange     :",profiler%mpi_halo_exchange,&
                        100.0*(profiler%mpi_halo_exchange/wall_clock)
                    WRITE(g_out,'(a23,2f16.4)')"Total                 :",kernel_total,&
                        100.0*(kernel_total/wall_clock)
                    WRITE(g_out,'(a23,2f16.4)')"The Rest              :",wall_clock-kernel_total,&
                        100.0*(wall_clock-kernel_total)/wall_clock
                ENDIF
            ENDIF

            EXIT

        END IF


    END DO


END SUBROUTINE hydro
