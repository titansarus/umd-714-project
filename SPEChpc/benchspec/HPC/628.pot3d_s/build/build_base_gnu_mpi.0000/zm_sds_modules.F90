!
! **********************************************************************
!
! Copyright 2018 Predictive Science Inc.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
! implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
! **********************************************************************
!
!#######################################################################
      module rp1d_def
!
!-----------------------------------------------------------------------
! ****** Define a structure to hold a REAL 1D pointer.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
      type :: rp1d
        real(r_typ), dimension(:), pointer :: f
      end type
!
      end module
!#######################################################################
      module sds_def
!
!-----------------------------------------------------------------------
! ****** Definition of the SDS data structure.
!-----------------------------------------------------------------------
!
      use number_types
      use rp1d_def
!
      implicit none
!
      integer, parameter, private :: mxdim=3
!
      type :: sds
        integer :: ndim
        integer, dimension(mxdim) :: dims
        logical :: scale
        logical :: hdf32
        type(rp1d), dimension(mxdim) :: scales
        real(r_typ), dimension(:,:,:), pointer :: f
      end type
!
      end module
!#######################################################################
      module rdhdf_1d_interface
      interface
        subroutine rdhdf_1d (fname,scale,nx,f,x,ierr)
        use number_types
        implicit none
        character(*) :: fname
        logical :: scale
        integer :: nx
        real(r_typ), dimension(:), pointer :: f
        real(r_typ), dimension(:), pointer :: x
        integer :: ierr
        end subroutine
      end interface
      end module
!#######################################################################
      module rdhdf_2d_interface
      interface
        subroutine rdhdf_2d (fname,scale,nx,ny,f,x,y,ierr)
        use number_types
        implicit none
        character(*) :: fname
        logical :: scale
        integer :: nx,ny
        real(r_typ), dimension(:,:), pointer :: f
        real(r_typ), dimension(:), pointer :: x,y
        integer :: ierr
        end subroutine
      end interface
      end module
!#######################################################################
      module rdhdf_3d_interface
      interface
        subroutine rdhdf_3d (fname,scale,nx,ny,nz,f,x,y,z,ierr)
        use number_types
        implicit none
        character(*) :: fname
        logical :: scale
        integer :: nx,ny,nz
        real(r_typ), dimension(:,:,:), pointer :: f
        real(r_typ), dimension(:), pointer :: x,y,z
        integer :: ierr
        end subroutine
      end interface
      end module
!#######################################################################
      module rdtxt_1d_interface
      interface
        subroutine rdtxt_1d (fname,scale,nx,f,x,ierr)
        use number_types
        implicit none
        character(*) :: fname
        logical :: scale
        integer :: nx
        real(r_typ), dimension(:), pointer :: f
        real(r_typ), dimension(:), pointer :: x
        integer :: ierr
        end subroutine
      end interface
      end module
!#######################################################################
      module rdtxt_2d_interface
      interface
        subroutine rdtxt_2d (fname,scale,nx,ny,f,x,y,ierr)
        use number_types
        implicit none
        character(*) :: fname
        logical :: scale
        integer :: nx,ny
        real(r_typ), dimension(:,:), pointer :: f
        real(r_typ), dimension(:), pointer :: x,y
        integer :: ierr
        end subroutine
      end interface
      end module
!#######################################################################
      module rdtxt_3d_interface
      interface
        subroutine rdtxt_3d (fname,scale,nx,ny,nz,f,x,y,z,ierr)
        use number_types
        implicit none
        character(*) :: fname
        logical :: scale
        integer :: nx,ny,nz
        real(r_typ), dimension(:,:,:), pointer :: f
        real(r_typ), dimension(:), pointer :: x,y,z
        integer :: ierr
        end subroutine
      end interface
      end module
