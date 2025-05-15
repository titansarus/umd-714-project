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

!>  @brief Fortran mpi buffer packing kernel
!>  @author Wayne Gaudin
!>  @details Packs/unpacks mpi send and receive buffers

MODULE pack_kernel_module

CONTAINS

  SUBROUTINE clover_pack_message_left(x_min,x_max,y_min,y_max,field,                &
                                      left_snd_buffer,                              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth,field_type,                             &
                                      buffer_offset)

    IMPLICIT NONE

    REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
    REAL(KIND=8) :: left_snd_buffer(:)

    INTEGER      :: CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA
    INTEGER      :: depth,field_type,x_min,x_max,y_min,y_max
    INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset

    ! Pack

    ! These array modifications still need to be added on, plus the donor data location changes as in update_halo
    IF(field_type.EQ.CELL_DATA) THEN
      x_inc=0
      y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
      x_inc=1
      y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
      x_inc=1
      y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
      x_inc=0
      y_inc=1
    ENDIF

#if defined(SPEC_OPENMP) 
    !$OMP PARALLEL DO PRIVATE(index)
#endif 
#ifdef SPEC_OPENMP_TARGET
!$OMP TARGET DATA MAP(to:left_snd_buffer, field)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(index)
#endif
#ifdef SPEC_OPENACC 
#  ifndef SPEC_ACCEL_AWARE_MPI
!$ACC DATA &
!$ACC PRESENT(left_snd_buffer, field)
#  endif
!$ACC KERNELS
!$ACC LOOP INDEPENDENT
#endif 
    DO k=y_min-depth,y_max+y_inc+depth
#ifdef SPEC_OPENACC 
!$ACC LOOP INDEPENDENT
#endif 
      DO j=1,depth
        index= buffer_offset + j+(k+depth-1)*depth
        left_snd_buffer(index)=field(x_min+x_inc-1+j,k)
      ENDDO
    ENDDO
#ifdef SPEC_OPENACC 
!$ACC END KERNELS
#  ifndef SPEC_ACCEL_AWARE_MPI
!$ACC UPDATE HOST(left_snd_buffer)
!$ACC END DATA
#  endif 
#endif 

#ifdef SPEC_OPENMP_TARGET
!$OMP TARGET UPDATE FROM(left_snd_buffer)
!$OMP END TARGET DATA
#endif

#if defined(SPEC_OPENMP) 
  !$OMP END PARALLEL DO
#endif 

  END SUBROUTINE clover_pack_message_left

  SUBROUTINE clover_unpack_message_left(x_min,x_max,y_min,y_max,field,                &
                                        left_rcv_buffer,                              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth,field_type,                             &
                                        buffer_offset)

    IMPLICIT NONE

    REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
    REAL(KIND=8) :: left_rcv_buffer(:)

    INTEGER      :: CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA
    INTEGER      :: depth,field_type,x_min,x_max,y_min,y_max
    INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset

    ! Unpack

    ! These array modifications still need to be added on, plus the donor data location changes as in update_halo
    IF(field_type.EQ.CELL_DATA) THEN
      x_inc=0
      y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
      x_inc=1
      y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
      x_inc=1
      y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
      x_inc=0
      y_inc=1
    ENDIF
#ifdef SPEC_OPENACC
#  ifndef SPEC_ACCEL_AWARE_MPI
!$ACC DATA &
!$ACC PRESENT (left_rcv_buffer, field)
!$ACC UPDATE DEVICE (left_rcv_buffer)
#  endif 
!$ACC KERNELS
!$ACC LOOP INDEPENDENT
#endif 

#ifdef SPEC_OPENMP_TARGET
!$OMP TARGET DATA MAP(to:left_rcv_buffer, field)
!$OMP TARGET UPDATE TO(left_rcv_buffer)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(index)
#endif

#if defined(SPEC_OPENMP) 
    !$OMP PARALLEL DO PRIVATE(index)
#endif 
    DO k=y_min-depth,y_max+y_inc+depth
#ifdef SPEC_OPENACC
!$ACC LOOP INDEPENDENT
#endif 
      DO j=1,depth
        index= buffer_offset + j+(k+depth-1)*depth
        field(x_min-j,k)=left_rcv_buffer(index)
      ENDDO
    ENDDO
#ifdef SPEC_OPENACC
!$ACC END KERNELS
#  ifndef SPEC_ACCEL_AWARE_MPI
!$ACC END DATA
#  endif 
#endif 
#ifdef SPEC_OPENMP_TARGET
!$OMP END TARGET DATA
#endif

#if defined(SPEC_OPENMP) 
  !$OMP END PARALLEL DO
#endif 

  END SUBROUTINE clover_unpack_message_left

  SUBROUTINE clover_pack_message_right(x_min,x_max,y_min,y_max,field,                &
                                       right_snd_buffer,                             &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth,field_type,                             &
                                       buffer_offset)

    IMPLICIT NONE

    REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
    REAL(KIND=8) :: right_snd_buffer(:)

    INTEGER      :: CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA
    INTEGER      :: depth,field_type,x_min,x_max,y_min,y_max
    INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset

    ! Pack

    ! These array modifications still need to be added on, plus the donor data location changes as in update_halo
    IF(field_type.EQ.CELL_DATA) THEN
      x_inc=0
      y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
      x_inc=1
      y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
      x_inc=1
      y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
      x_inc=0
      y_inc=1
    ENDIF
#ifdef SPEC_OPENACC 
#  ifndef SPEC_ACCEL_AWARE_MPI
!$ACC DATA &
!$ACC PRESENT(right_snd_buffer, field)
#  endif
!$ACC KERNELS
!$ACC LOOP INDEPENDENT
#endif  

#ifdef SPEC_OPENMP_TARGET
!$OMP TARGET DATA MAP(to:right_snd_buffer, field)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(index)
#endif

#if defined(SPEC_OPENMP) 
    !$OMP PARALLEL DO PRIVATE(index)
#endif 
    DO k=y_min-depth,y_max+y_inc+depth
#ifdef SPEC_OPENACC
!$ACC LOOP INDEPENDENT
#endif 
      DO j=1,depth
        index= buffer_offset + j+(k+depth-1)*depth
        right_snd_buffer(index)=field(x_max+1-j,k)
      ENDDO
    ENDDO
#ifdef SPEC_OPENACC
!$ACC END KERNELS
#  ifndef SPEC_ACCEL_AWARE_MPI
!$ACC UPDATE HOST(right_snd_buffer)
!$ACC END DATA
#  endif
#endif  

#ifdef SPEC_OPENMP_TARGET
!$OMP TARGET UPDATE FROM(right_snd_buffer)
!$OMP END TARGET DATA
#endif

#if defined(SPEC_OPENMP) 
  !$OMP END PARALLEL DO
#endif 

  END SUBROUTINE clover_pack_message_right

  SUBROUTINE clover_unpack_message_right(x_min,x_max,y_min,y_max,field,                &
                                         right_rcv_buffer,                             &
                                         CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                         depth,field_type,                             &
                                         buffer_offset)

    IMPLICIT NONE

    REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
    REAL(KIND=8) :: right_rcv_buffer(:)

    INTEGER      :: CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA
    INTEGER      :: depth,field_type,x_min,x_max,y_min,y_max
    INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset

    ! Unpack

    ! These array modifications still need to be added on, plus the donor data location changes as in update_halo
    IF(field_type.EQ.CELL_DATA) THEN
      x_inc=0
      y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
      x_inc=1
      y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
      x_inc=1
      y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
      x_inc=0
      y_inc=1
    ENDIF

#if defined(SPEC_OPENMP) 
    !$OMP PARALLEL DO PRIVATE(index)
#endif 

#ifdef SPEC_OPENMP_TARGET
!$OMP TARGET DATA MAP(to:right_rcv_buffer)
!$OMP TARGET UPDATE TO(right_rcv_buffer)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(index)
#endif

#ifdef SPEC_OPENACC
#  ifndef SPEC_ACCEL_AWARE_MPI
!$ACC DATA &
!$ACC PRESENT (right_rcv_buffer, field)
!$ACC UPDATE DEVICE (right_rcv_buffer)
#  endif
!$ACC KERNELS
!$ACC LOOP INDEPENDENT
#endif 
    DO k=y_min-depth,y_max+y_inc+depth
#ifdef SPEC_OPENACC
!$ACC LOOP INDEPENDENT
#endif 
      DO j=1,depth
        index= buffer_offset + j+(k+depth-1)*depth
        field(x_max+x_inc+j,k)=right_rcv_buffer(index)
      ENDDO
    ENDDO
#ifdef SPEC_OPENACC
!$ACC END KERNELS
#  ifndef SPEC_ACCEL_AWARE_MPI
!$ACC END DATA
#  endif
#endif 
#ifdef SPEC_OPENMP_TARGET
!$OMP END TARGET DATA
#endif 
#if defined(SPEC_OPENMP) 
  !$OMP END PARALLEL DO
#endif 

  END SUBROUTINE clover_unpack_message_right

  SUBROUTINE clover_pack_message_top(x_min,x_max,y_min,y_max,field,                &
                                     top_snd_buffer,                               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth,field_type,                             &
                                     buffer_offset)

    IMPLICIT NONE

    REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
    REAL(KIND=8) :: top_snd_buffer(:)

    INTEGER      :: CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA
    INTEGER      :: depth,field_type,x_min,x_max,y_min,y_max
    INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset

    ! Pack

    ! These array modifications still need to be added on, plus the donor data location changes as in update_halo
    IF(field_type.EQ.CELL_DATA) THEN
      x_inc=0
      y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
      x_inc=1
      y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
      x_inc=1
      y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
      x_inc=0
      y_inc=1
    ENDIF

#ifdef SPEC_OPENMP_TARGET
!$OMP TARGET DATA MAP(to:top_snd_buffer, field)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(index)
#endif

#ifdef SPEC_OPENACC
#  ifndef SPEC_ACCEL_AWARE_MPI
!$ACC DATA &
!$ACC PRESENT(top_snd_buffer, field)
#  endif
!$ACC KERNELS
!$ACC LOOP INDEPENDENT
#endif 
    DO k=1,depth
#if defined(SPEC_OPENMP) 
      !$OMP PARALLEL DO PRIVATE(index)
#endif 
#ifdef SPEC_OPENACC
!$ACC LOOP INDEPENDENT
#endif 
      DO j=x_min-depth,x_max+x_inc+depth
        index= buffer_offset + k+(j+depth-1)*depth
        top_snd_buffer(index)=field(j,y_max+1-k)
      ENDDO
#if defined(SPEC_OPENMP) 
    !$OMP END PARALLEL DO
#endif  
    ENDDO

#ifdef SPEC_OPENMP_TARGET
!$OMP TARGET UPDATE FROM(top_snd_buffer)
!$OMP END TARGET DATA
#endif
#ifdef SPEC_OPENACC
!$ACC END KERNELS
#  ifndef SPEC_ACCEL_AWARE_MPI
!$ACC UPDATE HOST(top_snd_buffer)
!$ACC END DATA
# endif
#endif 
  
  END SUBROUTINE clover_pack_message_top

  SUBROUTINE clover_unpack_message_top(x_min,x_max,y_min,y_max,field,                &
                                       top_rcv_buffer,                               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth,field_type,                             &
                                       buffer_offset)

    IMPLICIT NONE

    REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
    REAL(KIND=8) :: top_rcv_buffer(:)

    INTEGER      :: CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA
    INTEGER      :: depth,field_type,x_min,x_max,y_min,y_max
    INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset

    ! Unpack

    ! These array modifications still need to be added on, plus the donor data location changes as in update_halo
    IF(field_type.EQ.CELL_DATA) THEN
      x_inc=0
      y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
      x_inc=1
      y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
      x_inc=1
      y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
      x_inc=0
      y_inc=1
    ENDIF

#ifdef SPEC_OPENACC
#  ifndef SPEC_ACCEL_AWARE_MPI
!$ACC DATA &
!$ACC PRESENT (top_rcv_buffer, field)
!$ACC UPDATE DEVICE (top_rcv_buffer)
#  endif
!$ACC KERNELS
!$ACC LOOP INDEPENDENT
#endif  
#ifdef SPEC_OPENMP_TARGET
!$OMP TARGET DATA MAP(to:top_rcv_buffer, field)
!$OMP TARGET UPDATE TO(top_rcv_buffer)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(index)
#endif

    DO k=1,depth
#if defined(SPEC_OPENMP) 
      !$OMP PARALLEL DO PRIVATE(index)
#endif 
#ifdef SPEC_OPENACC
!$ACC LOOP INDEPENDENT
#endif 
      DO j=x_min-depth,x_max+x_inc+depth
        index= buffer_offset + k+(j+depth-1)*depth
        !index= buffer_offset + j + depth+(k-1)*(x_max+x_inc+(2*depth))
        field(j,y_max+y_inc+k)=top_rcv_buffer(index)
      ENDDO
#if defined(SPEC_OPENMP) 
    !$OMP END PARALLEL DO
#endif 
    ENDDO
#ifdef SPEC_OPENACC
!$ACC END KERNELS
#  ifndef SPEC_ACCEL_AWARE_MPI
!$ACC END DATA
#  endif
#endif 
#ifdef SPEC_OPENMP_TARGET
!$OMP END TARGET DATA
#endif

  END SUBROUTINE clover_unpack_message_top

  SUBROUTINE clover_pack_message_bottom(x_min,x_max,y_min,y_max,field,                &
                                        bottom_snd_buffer,                            &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth,field_type,                             &
                                        buffer_offset)

    IMPLICIT NONE

    REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
    REAL(KIND=8) :: bottom_snd_buffer(:)

    INTEGER      :: CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA
    INTEGER      :: depth,field_type,x_min,x_max,y_min,y_max
    INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset

    ! Pack

    ! These array modifications still need to be added on, plus the donor data location changes as in update_halo
    IF(field_type.EQ.CELL_DATA) THEN
      x_inc=0
      y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
      x_inc=1
      y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
      x_inc=1
      y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
      x_inc=0
      y_inc=1
    ENDIF

#ifdef SPEC_OPENACC
#  ifndef SPEC_ACCEL_AWARE_MPI
!$ACC DATA &
!$ACC PRESENT(bottom_snd_buffer, field)
#  endif
!$ACC KERNELS
!$ACC LOOP INDEPENDENT
#endif 

#ifdef SPEC_OPENMP_TARGET
!$OMP TARGET DATA MAP(to:bottom_snd_buffer, field)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(index)
#endif
    DO k=1,depth
#if defined(SPEC_OPENMP) 
      !$OMP PARALLEL DO PRIVATE(index)
#endif 
#ifdef SPEC_OPENACC
!$ACC LOOP INDEPENDENT
#endif 
      DO j=x_min-depth,x_max+x_inc+depth
        index= buffer_offset + k+(j+depth-1)*depth
        !index= buffer_offset + j+depth+(k-1)*(x_max+x_inc+(2*depth))
        bottom_snd_buffer(index)=field(j,y_min+y_inc-1+k)
      ENDDO
#if defined(SPEC_OPENMP) 
    !$OMP END PARALLEL DO
#endif 
    ENDDO
#ifdef SPEC_OPENACC
!$ACC END KERNELS
#  ifndef SPEC_ACCEL_AWARE_MPI
!$ACC UPDATE HOST(bottom_snd_buffer)
!$ACC END DATA
#  endif
#endif 

#ifdef SPEC_OPENMP_TARGET
!$OMP TARGET UPDATE FROM(bottom_snd_buffer)
!$OMP END TARGET DATA
#endif

  END SUBROUTINE clover_pack_message_bottom

  SUBROUTINE clover_unpack_message_bottom(x_min,x_max,y_min,y_max,field,                &
                                          bottom_rcv_buffer,                            &
                                          CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                          depth,field_type,                             &
                                          buffer_offset)

    IMPLICIT NONE

    REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
    REAL(KIND=8) :: bottom_rcv_buffer(:)

    INTEGER      :: CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA
    INTEGER      :: depth,field_type,x_min,x_max,y_min,y_max
    INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset

    ! Unpack

    ! These array modifications still need to be added on, plus the donor data location changes as in update_halo
    IF(field_type.EQ.CELL_DATA) THEN
      x_inc=0
      y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
      x_inc=1
      y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
      x_inc=1
      y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
      x_inc=0
      y_inc=1
    ENDIF

#if defined(SPEC_OPENMP) 
    !$OMP PARALLEL
#endif 
#ifdef SPEC_OPENACC
#  ifndef SPEC_ACCEL_AWARE_MPI
!$ACC DATA &
!$ACC PRESENT (bottom_rcv_buffer, field)
!$ACC UPDATE DEVICE (bottom_rcv_buffer)
#  endif
!$ACC KERNELS
!$ACC LOOP INDEPENDENT
#endif 
#ifdef SPEC_OPENMP_TARGET
!$OMP TARGET DATA MAP(to:bottom_rcv_buffer, field)
!$OMP TARGET UPDATE TO(bottom_rcv_buffer)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(index)
#endif

    DO k=1,depth
#if defined(SPEC_OPENMP)
      !$OMP DO PRIVATE(index)
#endif 
#ifdef SPEC_OPENACC
!$ACC LOOP INDEPENDENT
#endif 
      DO j=x_min-depth,x_max+x_inc+depth
        index= buffer_offset + k+(j+depth-1)*depth
        !index= buffer_offset + j+depth+(k-1)*(x_max+x_inc+(2*depth))
        field(j,y_min-k)=bottom_rcv_buffer(index)
      ENDDO
#if defined(SPEC_OPENMP) 
    !$OMP END DO
#endif 
    ENDDO
#ifdef SPEC_OPENACC
!$ACC END KERNELS
#  ifndef SPEC_ACCEL_AWARE_MPI
!$ACC END DATA
#  endif
#endif 
#ifdef SPEC_OPENMP_TARGET
!$OMP END TARGET DATA
#endif 
#if defined(SPEC_OPENMP) 
  !$OMP END PARALLEL
#endif 

  END SUBROUTINE clover_unpack_message_bottom

END MODULE pack_kernel_module
