!
!-----------------------------------------------------------------------
!
! ****** Source to build the SDS library.
! ****** These routines are used by Zoran Mikic's tools.
!
!-----------------------------------------------------------------------
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
!        07/29/2003, ZM, Version 1.00:
!
!         - Original version of the SDS library.
!           This library was put together to facilitate the
!           development of ZM's tools.
!           It includes the new read/write routines for
!           scientific data sets (both text and HDF format).
!           The code was cleaned up to use standard FORTRAN90.
!
!        02/20/2004, ZM, Version 1.01:
!
!         - Added the ability to specify the format in writing
!           floating point numbers in routine WRFP.  This is used
!           in writing SDS text files using routine WRTXT.
!           For 32-bit data, WRTXT specifies 7 digits of precision,
!           whereas for 64-bit data, 15 digits of precision are
!           used.
!
!        04/02/2005, ZM, Version 1.02:
!
!         - Added a call to DFSDclear() in routine WRHDF.  This
!           initializes the SDS interface to the default state
!           for each new file.  This is needed to prevent settings
!           from previous files from interfering with each other.
!
!        06/16/2006, ZM, Version 1.03:
!
!         - Fixed some pointer allocation and manipulation issues in
!           the HDF read and write routines to make them obey
!           the FORTRAN 90 standard.  They were misbehaving on
!           the IFORT compiler.
!
!        02/24/2009, ZM, Version 1.04:
!
!         - Made a small change to the way an SDS is deallocated.
!         - Added a routine to initialize an SDS.  This is useful
!           when deallocating SDS structures.
!
!        08/30/2016, RC, Version 2.00:
!
!         - Added ability to read and write hdf5 files.
!           Now, rdhdf and wrhdf will read or write an hdf5 file
!           if the given fname ends in ".h5".
!           The library now needs to be linked to the hdf5 libraries.
!         - Modified rdhdf to be compatible with hdf4 files made using
!           the SD API instead of the DFSD API.
!
!        09/05/2016, RC, Version 2.01:
!
!         - Fixed problem with hdf5 writes when using 32-bit data.
!
!        05/22/2017, RC, Version 2.02:
!
!         - Fixed problem with 1D and 2D hdf5 writes when
!           the s%dims() are not set to 1 for the unused dimensions.
!
!        05/03/2019, RC, Version 2.03:
!
!         - Bug fix for HDF5 reads.
!-----------------------------------------------------------------------
!
!#######################################################################
      module sdslib_ident
!
      character(*), parameter :: cname='SDSLIB'
      character(*), parameter :: cvers='2.03'
      character(*), parameter :: cdate='05/03/2019'
!
      end module
!#######################################################################
      module assign_ptr_1d_interface
      interface
        subroutine assign_ptr_1d (from,to)
        use number_types
        implicit none
        real(r_typ), dimension(:), target :: from
        real(r_typ), dimension(:), pointer :: to
        end subroutine
      end interface
      end module
!#######################################################################
      module assign_ptr_3d_interface
      interface
        subroutine assign_ptr_3d (from,to)
        use number_types
        implicit none
        real(r_typ), dimension(:,:,:), target :: from
        real(r_typ), dimension(:,:,:), pointer :: to
        end subroutine
      end interface
      end module
!#######################################################################
      subroutine assign_ptr_1d (from,to)
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(:), target :: from
      real(r_typ), dimension(:), pointer :: to
!
!-----------------------------------------------------------------------
!
      to=>from
!
      return
      end subroutine
!#######################################################################
      subroutine assign_ptr_3d (from,to)
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(:,:,:), target :: from
      real(r_typ), dimension(:,:,:), pointer :: to
!
!-----------------------------------------------------------------------
!
      to=>from
!
      return
      end subroutine
!#######################################################################
      subroutine init_sds_pointer_status (s)
!
!-----------------------------------------------------------------------
!
! ****** Disassociate all the pointers in the SDS in structure S.
!
!-----------------------------------------------------------------------
!
! ****** This is useful when subsequently querying the association
! ****** status of these pointers (e.g., when deallocating storage).
!
!-----------------------------------------------------------------------
!
      use sds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(sds) :: s
!
!-----------------------------------------------------------------------
!
      nullify (s%f)
!
      nullify (s%scales(1)%f)
      nullify (s%scales(2)%f)
      nullify (s%scales(3)%f)
!
      return
      end subroutine
!#######################################################################
      subroutine deallocate_sds (s)
!
!-----------------------------------------------------------------------
!
! ****** Deallocate the memory used by the SDS in structure S.
!
!-----------------------------------------------------------------------
!
      use sds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(sds) :: s
!
!-----------------------------------------------------------------------
!
      if (associated(s%f)) deallocate (s%f)
!
      if (associated(s%scales(1)%f)) deallocate (s%scales(1)%f)
      if (associated(s%scales(2)%f)) deallocate (s%scales(2)%f)
      if (associated(s%scales(3)%f)) deallocate (s%scales(3)%f)
!
      return
      end subroutine
!#######################################################################
      subroutine rdhdf_1d (fname,scale,nx,f,x,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read a 1D scientific data set from an HDF file.
!
!-----------------------------------------------------------------------
!
! ****** This routine calls routine RDHDF to read the file.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      logical :: scale
      integer :: nx
      real(r_typ), dimension(:), pointer :: f
      real(r_typ), dimension(:), pointer :: x
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ierr
!
!-----------------------------------------------------------------------
!
! ****** Declaration for the SDS structure.
!
      type(sds) :: s
!
!-----------------------------------------------------------------------
!
! ****** Read the data set.
!
      call rdhdf (fname,s,ierr)
!
      if (ierr.ne.0) return
!
! ****** Check that this is a 1D data set.
!
      if (s%ndim.ne.1) then
        write (*,*)
        write (*,*) '### ERROR in RDHDF_1D:'
        write (*,*) '### The HDF file does not contain a 1D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
!
! ****** Set the output arguments.
!
      nx=s%dims(1)
      scale=s%scale
      x=>s%scales(1)%f
!
      allocate (f(nx))
      f=s%f(:,1,1)
      deallocate (s%f)
!
      return
      end
!#######################################################################
      subroutine rdhdf_2d (fname,scale,nx,ny,f,x,y,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read a 2D scientific data set from an HDF file.
!
!-----------------------------------------------------------------------
!
! ****** This routine calls routine RDHDF to read the file.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      logical :: scale
      integer :: nx,ny
      real(r_typ), dimension(:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: x,y
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ny,ierr
!
!-----------------------------------------------------------------------
!
! ****** Declaration for the SDS structure.
!
      type(sds) :: s
!
!-----------------------------------------------------------------------
!
! ****** Read the data set.
!
      call rdhdf (fname,s,ierr)
!
      if (ierr.ne.0) return
!
! ****** Check that this is a 2D data set.
!
      if (s%ndim.ne.2) then
        write (*,*)
        write (*,*) '### ERROR in RDHDF_2D:'
        write (*,*) '### The HDF file does not contain a 2D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
!
! ****** Set the output arguments.
!
      nx=s%dims(1)
      ny=s%dims(2)
      scale=s%scale
      x=>s%scales(1)%f
      y=>s%scales(2)%f
!
      allocate (f(nx,ny))
      f=s%f(:,:,1)
      deallocate (s%f)
!
      return
      end
!#######################################################################
      subroutine rdhdf_3d (fname,scale,nx,ny,nz,f,x,y,z,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read a 3D scientific data set from an HDF file.
!
!-----------------------------------------------------------------------
!
! ****** This routine calls routine RDHDF to read the file.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      logical :: scale
      integer :: nx,ny,nz
      real(r_typ), dimension(:,:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: x,y,z
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ny,nz,ierr
!
!-----------------------------------------------------------------------
!
! ****** Declaration for the SDS structure.
!
      type(sds) :: s
!
!-----------------------------------------------------------------------
!
! ****** Read the data set.
!
      call rdhdf (fname,s,ierr)
!
      if (ierr.ne.0) return
!
! ****** Check that this is a 3D data set.
!
      if (s%ndim.ne.3) then
        write (*,*)
        write (*,*) '### ERROR in RDHDF_3D:'
        write (*,*) '### The HDF file does not contain a 3D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
!
! ****** Set the output arguments.
!
      nx=s%dims(1)
      ny=s%dims(2)
      nz=s%dims(3)
      scale=s%scale
      x=>s%scales(1)%f
      y=>s%scales(2)%f
      z=>s%scales(3)%f
      f=>s%f
!
      return
      end
!#######################################################################
      subroutine wrhdf_1d (fname,scale,nx,f,x,hdf32,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Write a 1D scientific data set to an HDF file.
!
!-----------------------------------------------------------------------
!
! ****** This routine calls routine WRHDF to write the file.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      logical :: scale
      integer :: nx
      real(r_typ), dimension(nx,1,1) :: f
      real(r_typ), dimension(nx) :: x
      logical :: hdf32
      integer :: ierr
      intent(in) :: fname,scale,nx,f,x,hdf32
      intent(out) :: ierr
!
!-----------------------------------------------------------------------
!
! ****** Declaration for the SDS structure.
!
      type(sds) :: s
!
!-----------------------------------------------------------------------
!
! ****** Set the structure components.
!
      s%ndim=1
      s%dims(1)=nx
      s%dims(2)=1
      s%dims(3)=1
      s%scale=scale
      s%hdf32=hdf32
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
      else
        nullify (s%scales(1)%f)
      end if
      nullify (s%scales(2)%f)
      nullify (s%scales(3)%f)
      call assign_ptr_3d (f,s%f)
!
! ****** Write the data set.
!
      call wrhdf (fname,s,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRHDF_1D:'
        write (*,*) '### Could not write the 1D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
!
      return
      end
!#######################################################################
      subroutine wrhdf_2d (fname,scale,nx,ny,f,x,y,hdf32,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Write a 2D scientific data set to an HDF file.
!
!-----------------------------------------------------------------------
!
! ****** This routine calls routine WRHDF to write the file.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      logical :: scale
      integer :: nx,ny
      real(r_typ), dimension(nx,ny,1) :: f
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      logical :: hdf32
      integer :: ierr
      intent(in) :: fname,scale,nx,ny,f,x,y,hdf32
      intent(out) :: ierr
!
!-----------------------------------------------------------------------
!
! ****** Declaration for the SDS structure.
!
      type(sds) :: s
!
!-----------------------------------------------------------------------
!
! ****** Set the structure components.
!
      s%ndim=2
      s%dims(1)=nx
      s%dims(2)=ny
      s%dims(3)=1
      s%scale=scale
      s%hdf32=hdf32
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
        call assign_ptr_1d (y,s%scales(2)%f)
      else
        nullify (s%scales(1)%f)
        nullify (s%scales(2)%f)
      end if
      nullify (s%scales(3)%f)
      call assign_ptr_3d (f,s%f)
!
! ****** Write the data set.
!
      call wrhdf (fname,s,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRHDF_2D:'
        write (*,*) '### Could not write the 2D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
!
      return
      end
!#######################################################################
      subroutine wrhdf_3d (fname,scale,nx,ny,nz,f,x,y,z,hdf32,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Write a 3D scientific data set to an HDF file.
!
!-----------------------------------------------------------------------
!
! ****** This routine calls routine WRHDF to write the file.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      logical :: scale
      integer :: nx,ny,nz
      real(r_typ), dimension(nx,ny,nz) :: f
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nz) :: z
      logical :: hdf32
      integer :: ierr
      intent(in) :: fname,scale,nx,ny,nz,f,x,y,z,hdf32
      intent(out) :: ierr
!
!-----------------------------------------------------------------------
!
! ****** Declaration for the SDS structure.
!
      type(sds) :: s
!
!-----------------------------------------------------------------------
!
! ****** Set the structure components.
!
      s%ndim=3
      s%dims(1)=nx
      s%dims(2)=ny
      s%dims(3)=nz
      s%scale=scale
      s%hdf32=hdf32
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
        call assign_ptr_1d (y,s%scales(2)%f)
        call assign_ptr_1d (z,s%scales(3)%f)
      else
        nullify (s%scales(1)%f)
        nullify (s%scales(2)%f)
        nullify (s%scales(3)%f)
      end if
      call assign_ptr_3d (f,s%f)
!
! ****** Write the data set.
!
      call wrhdf (fname,s,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRHDF_3D:'
        write (*,*) '### Could not write the 3D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
!
      return
      end
!#######################################################################
      subroutine rdsds (fmt,fname,s,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read a scientific data set from file FNAME into
! ****** SDS structure S.
!
! ****** Use routine RDTXT or RDHDF, depending on the format
! ****** specified by FMT.
!
!-----------------------------------------------------------------------
!
      use sds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fmt
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fmt,fname
      intent(out) :: s,ierr
!
!-----------------------------------------------------------------------
!
      if (fmt.eq.'text') then
        call rdtxt (fname,s,ierr)
      else if (fmt.eq.'hdf'.or.fmt.eq.'h5') then
        call rdhdf (fname,s,ierr)
      else
        write (*,*)
        write (*,*) '### ERROR in RDSDS:'
        write (*,*) '### Invalid file format specified.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) 'Format = ',fmt
        ierr=5
        return
      end if
!
      return
      end
!#######################################################################
      subroutine wrsds (fmt,fname,s,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Write a scientific data set from SDS structure S to
! ****** file FNAME.
!
! ****** Use routine WRTXT or WRHDF, depending on the format
! ****** specified by FMT.
!
!-----------------------------------------------------------------------
!
      use sds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fmt
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fmt,fname,s
      intent(out) :: ierr
!
!-----------------------------------------------------------------------
!
      if (fmt.eq.'text') then
        call wrtxt (fname,s,ierr)
      else if (fmt.eq.'hdf'.or.fmt.eq.'h5') then
        call wrhdf (fname,s,ierr)
      else
        write (*,*)
        write (*,*) '### ERROR in WRSDS:'
        write (*,*) '### Invalid file format specified.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) 'Format = ',fmt
        ierr=5
        return
      end if
!
      return
      end
!#######################################################################
      subroutine rdhdf (fname,s,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read a 1D, 2D, or 3D scientific data set from an HDF file.
! ****** This routine uses the new SD API instead of the
! ****** outdated DFSD API.
!
!-----------------------------------------------------------------------
!
! ****** This routine allocates the required memory and returns
! ****** pointers to the data and scale arrays.
!
!-----------------------------------------------------------------------
!
! ****** Input arguments:
!
!          FNAME   : [character(*)]
!                    HDF data file name to read from.
!
! ****** Output arguments:
!
!          S       : [structure of type SDS]
!                    A structure that holds the field, its
!                    dimensions, and the scales, with the
!                    components described below.
!
!          IERR    : [integer]
!                    IERR=0 is returned if the data set was read
!                    successfully.  Otherwise, IERR is set to a
!                    nonzero value.
!
! ****** Components of structure S:
!
!          NDIM    : [integer]
!                    Number of dimensions found in the data set.
!
!          DIMS    : [integer, dimension(3)]
!                    Number of points in the data set dimensions.
!                    For a 1D data set, DIMS(2)=DIMS(3)=1.
!                    For a 2D data set, DIMS(3)=1.
!
!          SCALE   : [logical]
!                    Flag to indicate the presence of scales (axes)
!                    in the data set.  SCALE=.false. means that scales
!                    were not found; SCALE=.true. means that scales
!                    were found.
!
!          HDF32   : [logical]
!                    Flag to indicate the precision of the data set
!                    read in.  HDF32=.true. means that the data is
!                    32-bit; HDF32=.false. means that the data is
!                    64-bit.
!
!          SCALES  : [structure of type RP1D, dimension(3)]
!                    This array holds the pointers to the scales
!                    when SCALE=.true., and is undefined otherwise.
!
!          F       : [real, pointer to a rank-3 array]
!                    This array holds the data set values.
!
! ****** The storage for the arrays pointed to by F, and the
! ****** scales (if present) in structure SCALES, is allocated by
! ****** this routine.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      character, dimension(256) :: sds_name, dim_name
      type(sds) :: s
      integer :: ierr,i
      intent(in) :: fname
      intent(out) :: s,ierr
!
!-----------------------------------------------------------------------
!
      ierr=0
!
!-----------------------------------------------------------------------
!
! ****** Read hdf5 file if fname ends in '.h5'.
!
      i=index(fname,'.h');
      if (fname(i+1:i+2).eq.'h5') then
        call rdh5 (fname,s,ierr)
        return
      else
        print*,"HDF4 has been disabled"
        ierr=-1
      end if
!
      return
      end
!#######################################################################
      subroutine rdh5 (fname,s,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read a 1D, 2D, or 3D scientific data set from an HDF5 file.
! ****** The HDF5 file is currently assumed to contain only one
! ****** dataset (1D,2d,3D), with or without scales, in group "/",
! ****** and has no other data members.
!
!-----------------------------------------------------------------------
!
! ****** This routine allocates the required memory and returns
! ****** pointers to the data and scale arrays.
!
!-----------------------------------------------------------------------
!
! ****** Input arguments:
!
!          FNAME   : [character(*)]
!                    HDF5 data file name to read from.
!
! ****** Output arguments:
!
!          S       : [structure of type SDS]
!                    A structure that holds the field, its
!                    dimensions, and the scales, with the
!                    components described below.
!
!          IERR    : [integer]
!                    IERR=0 is returned if the data set was read
!                    successfully.  Otherwise, IERR is set to a
!                    nonzero value.
!
! ****** Components of structure S:
!
!          NDIM    : [integer]
!                    Number of dimensions found in the data set.
!
!          DIMS    : [integer, dimension(3)]
!                    Number of points in the data set dimensions.
!                    For a 1D data set, DIMS(2)=DIMS(3)=1.
!                    For a 2D data set, DIMS(3)=1.
!
!          SCALE   : [logical]
!                    Flag to indicate the presence of scales (axes)
!                    in the data set.  SCALE=.false. means that scales
!                    were not found; SCALE=.true. means that scales
!                    were found.
!
!          HDF32   : [logical]
!                    Flag to indicate the precision of the data set
!                    read in.  HDF32=.true. means that the data is
!                    32-bit; HDF32=.false. means that the data is
!                    64-bit.
!
!          SCALES  : [structure of type RP1D, dimension(3)]
!                    This array holds the pointers to the scales
!                    when SCALE=.true., and is undefined otherwise.
!
!          F       : [real, pointer to a rank-3 array]
!                    This array holds the data set values.
!
! ****** The storage for the arrays pointed to by F, and the
! ****** scales (if present) in structure SCALES, is allocated by
! ****** this routine.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
      use hdf5
      use h5ds
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(sds) :: s
      character(*) :: fname
!
!-----------------------------------------------------------------------
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      integer :: i,obj_type,n_members
!
      integer(HID_T) :: file_id       ! File identifier
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(HID_T) :: dspace_id     ! Dataspace identifier
      integer(HID_T) :: dim_id        ! Dimension identifiers
      integer(HID_T) :: datatype_id   ! Datatype identifiers
!
      integer(SIZE_T) :: prec
!
      integer(HSIZE_T),dimension(:), allocatable :: s_dims,maxpts
      integer(HSIZE_T),dimension(1) :: s_dims_i
!
      real(KIND_REAL_4), dimension(:,:,:), allocatable :: f4
      real(KIND_REAL_4), dimension(:),     allocatable :: f4dim
      real(KIND_REAL_8), dimension(:,:,:), allocatable :: f8
      real(KIND_REAL_8), dimension(:),     allocatable :: f8dim
!
      character(512) :: obj_name
      character(4), parameter :: cname='RDH5'
!
      logical :: is_scale
!
!-----------------------------------------------------------------------
!
! ****** Initialize dimension count and arrays.
!
      s%ndim=0
      s%dims(:)=1
!
! ****** Initialize hdf5 interface.
!
      call h5open_f (ierr)
!
! ****** Open hdf5 file.
!
      call h5Fopen_f (trim(fname),H5F_ACC_RDONLY_F,file_id,ierr)
!
! ****** Get information about the hdf5 file.
!
      call h5Gn_members_f (file_id,"/",n_members,ierr)
!
! ****** Make sure there is (at maximum) one 3D dataset with scales.
!
      if (n_members.eq.0.or.n_members.gt.4) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Input file contains too few/many datasets.'
        write (*,*) 'File name: ',trim(fname)
        return
      endif
!
! ****** Assume the Dataset is in index 0 and get its name.
!
      call h5Gget_obj_info_idx_f (file_id,"/",0,obj_name,obj_type,ierr)
!
! ****** Open Dataset.
!
      call h5Dopen_f (file_id,trim(obj_name),dset_id,ierr)
!
! ****** Make sure the Dataset is not a scale.
!
      call h5DSis_scale_f(dset_id,is_scale,ierr)
      if (is_scale) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Input file Dataset at index 0 is a scale.'
        write (*,*) 'File name: ',trim(fname)
        return
      endif
!
! ****** Get dimensions (need s_dims array for format requirements).
!
      call h5Dget_space_f (dset_id,dspace_id,ierr)
      call h5Sget_simple_extent_ndims_f (dspace_id,s%ndim,ierr)
!
      allocate(s_dims(s%ndim))
!  
      allocate(maxpts(s%ndim))
      call h5Sget_simple_extent_dims_f (dspace_id,s_dims,maxpts,ierr)
      deallocate(maxpts)
!
      s%dims(1:s%ndim)=s_dims(1:s%ndim)
!
! ****** Get the floating-point precision of the data and set flag.
!
      call h5Dget_type_f (dset_id,datatype_id,ierr)
      call h5Tget_precision_f (datatype_id,prec,ierr)
!
      if (prec.eq.32) then
        s%hdf32=.true.
      elseif (prec.eq.64) then
        s%hdf32=.false.
      end if
!
! ****** Allocate the memory for the Dataset array in s.
!
      allocate (s%f(s%dims(1),s%dims(2),s%dims(3)))
!
! ****** Need to read the file in its own datatype, and then convert
! ****** to datatype of s%f.
!
      if (s%hdf32) then
        allocate (f4(s%dims(1),s%dims(2),s%dims(3)))
        call h5Dread_f (dset_id,datatype_id,f4,s_dims,ierr)
        s%f(:,:,:)=f4(:,:,:)
        deallocate (f4)
      else
        allocate (f8(s%dims(1),s%dims(2),s%dims(3)))
        call h5Dread_f (dset_id,datatype_id,f8,s_dims,ierr)
        s%f(:,:,:)=f8(:,:,:)
        deallocate (f8)
      end if
!
      deallocate(s_dims)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in RDH5:'
        write (*,*) '### Error while reading the dataset.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) '[Error return (from H5DREAD_F) = ',ierr,']'
        ierr=4
        return
      end if
!
! ****** Close the hdf5 type descriptor.
!
      call h5Tclose_f (datatype_id,ierr)
!
! ****** Check if there might be scales present, if so, read them.
!
      if (n_members.gt.1) then
!
! ***** First check that the number of scale datasets match the # dim.
!
        if (n_members-1.ne.s%ndim) then
          write (*,*)
          write (*,*) '### ERROR in RDH5:'
          write (*,*) '### # scales does not match # dims.'
          write (*,*) 'File name: ',trim(fname)
          return
        end if
!
        s%scale=.true.
!
! ****** Loop through scales, make sure each is a scale, and read them.
!
        do i=1,n_members-1
!
! ****** Get the name of scale dataset.
!
          call h5Gget_obj_info_idx_f (file_id,"/",i,obj_name,obj_type,ierr)
!
! ****** Open scale dataset.
!
          call h5Dopen_f (file_id,trim(obj_name),dim_id,ierr)
!
! ****** Make sure the scale is a scale.
!
          call h5DSis_scale_f (dim_id,is_scale,ierr)
          if (.not.is_scale) then
            write (*,*)
            write (*,*) '### ERROR in RDH5:'
            write (*,*) '### Scale is not a scale.'
            write (*,*) 'File name: ',trim(fname)
            return
          end if
!
! ****** Get dimension of scale.
!
          s_dims_i=s%dims(i)
!
! ****** Allocate scale.
!
          allocate (s%scales(i)%f(s_dims_i(1)))
!
! ****** Get the floating-point precision of the scale.
!
          call h5Dget_type_f (dim_id,datatype_id,ierr)
          call h5Tget_precision_f (datatype_id,prec,ierr)
!
! ****** Read in the scale data.
!
          if (s%hdf32) then
            allocate (f4dim(s_dims_i(1)))
            call h5Dread_f (dim_id,datatype_id,f4dim,s_dims_i,ierr)
            s%scales(i)%f(:)=f4dim(:)
            deallocate (f4dim)
          else
            allocate (f8dim(s_dims_i(1)))
            call h5Dread_f (dim_id,datatype_id,f8dim,s_dims_i,ierr)
            s%scales(i)%f(:)=f8dim(:)
            deallocate (f8dim)
          end if
!
! ****** Close the scale dataset.
!
          call h5Dclose_f (dim_id,ierr)
!
        enddo
!
! ****** Allocate dummy scales (of length 1) for empty dimensions.
!
        do i=s%ndim+1,3
          allocate (s%scales(i)%f(1))
        enddo
      else
!
! ****** If scales are not present, allocate dummy
! ****** scales (of length 1) so that the pointers to the scales
! ****** are valid.
!
        s%scale=.false.
!
        allocate (s%scales(1)%f(1))
        allocate (s%scales(2)%f(1))
        allocate (s%scales(3)%f(1))
      end if
!
! ****** Close the dataset.
!
      call h5Dclose_f (dset_id,ierr)
!
! ****** Close the file.
!
      call h5Fclose_f (file_id,ierr)
!
! ****** Close FORTRAN interface.
!
      call h5close_f (ierr)
!
      return
      end subroutine
!#######################################################################
      subroutine wrhdf (fname,s,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Write a 1D, 2D, or 3D scientific data set to an HDF file.
!
!-----------------------------------------------------------------------
!
! ****** Input arguments:
!
!          FNAME   : [character(*)]
!                    HDF data file name to write to.
!
!          S       : [structure of type SDS]
!                    A structure that holds the field, its
!                    dimensions, and the scales, with the
!                    components described below.
!
! ****** Output arguments:
!
!          IERR    : [integer]
!                    IERR=0 is returned if the data set was written
!                    successfully.  Otherwise, IERR is set to a
!                    nonzero value.
!
! ****** Components of structure S:
!
!          NDIM    : [integer]
!                    Number of dimensions in the data set.
!
!          DIMS    : [integer, dimension(3)]
!                    Number of points in the data set dimensions.
!                    Only DIMS(1 .. NDIM) are referenced.
!
!          SCALE   : [logical]
!                    Flag to indicate the presence of scales (axes)
!                    in the data set.  SCALE=.false. means that scales
!                    are not being supplied; SCALE=.true. means that
!                    scales are being supplied.
!
!          HDF32   : [logical]
!                    Flag to specify the precision of the data to
!                    be written to the file.  Set HDF32=.true. to
!                    write 32-bit data, and HDF32=.false. to write
!                    64-bit data.
!
!          SCALES  : [structure of type RP1D, dimension(3)]
!                    This array holds the pointers to the scales
!                    when SCALE=.true., and is not referenced
!                    otherwise.
!
!          F       : [real, pointer to a rank-3 array]
!                    This array holds the data set values.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fname,s
      intent(out) :: ierr
!
!-----------------------------------------------------------------------
!
      integer :: iret,i,n
!
!-----------------------------------------------------------------------
!
      ierr=0
!
! ****** Check the number of dimensions.
!
      if (s%ndim.le.0.or.s%ndim.gt.3) then
        write (*,*)
        write (*,*) '### ERROR in WRHDF:'
        write (*,*) '### Could not write the SDS data.'
        write (*,*) 'Invalid number of dimensions.'
        write (*,*) 'Number of dimensions = ',s%ndim
        write (*,*) 'File name: ',trim(fname)
        ierr=1
        return
      end if
!
! ****** Write hdf5 file if fname ends in '.h5'.
!
      i=index(fname,'.h')
      if (fname(i+1:i+2).eq.'h5') then
        call wrh5 (fname,s,ierr)
      else
        print*,"HDF4 disabled."
        ierr=-1
      end if
!
      return
      end
!#######################################################################
      subroutine wrh5 (fname,s,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Write a 1D, 2D, or 3D scientific data set to an HDF5 file.
!
!-----------------------------------------------------------------------
!
! ****** Input arguments:
!
!          FNAME   : [character(*)]
!                    HDF data file name to write to.
!
!          S       : [structure of type SDS]
!                    A structure that holds the field, its
!                    dimensions, and the scales, with the
!                    components described below.
!
! ****** Output arguments:
!
!          IERR    : [integer]
!                    IERR=0 is returned if the data set was written
!                    successfully.  Otherwise, IERR is set to a
!                    nonzero value.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
      use hdf5
      use h5ds
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fname,s
      intent(out) :: ierr
!
!-----------------------------------------------------------------------
!
      character(8) ::   dimname
      integer :: i
      integer(HID_T) :: file_id       ! File identifier
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(HID_T) :: dspace_id,dspacedim_id   ! Dataspace identifiers
      integer(HID_T) :: dim_id        ! Dimension identifiers
      integer(HSIZE_T),dimension(3) :: s_dims
      integer(HSIZE_T),dimension(1) :: s_dims_i
!
      real(KIND_REAL_4), dimension(:,:,:), allocatable :: f4
      real(KIND_REAL_4), dimension(:),     allocatable :: f4dim
      real(KIND_REAL_8), dimension(:,:,:), allocatable :: f8
      real(KIND_REAL_8), dimension(:),     allocatable :: f8dim
!
!-----------------------------------------------------------------------
!
! ****** HDF5 calls are picky about the integer format for the dims
! ****** so the s%dims need to be converted to HSIZE_T integers.
!
! ****** Also, sometimes calls to wrhdf() for 1D and 2D datasets
! ****** do not have the unused dims(i) set to 1 (unset).
! ****** To avoid needing a function call to implicitly reshape
! ****** f(n), set the dims here.
!
      do i=1,3
         if (i.le.s%ndim) then
           s_dims(i)=s%dims(i)
         else
           s_dims(i)=1
         endif
      end do
!
! ****** Initialize hdf5 interface.
!
      call h5open_f (ierr)
!
! ****** Create the file.
!
      call h5Fcreate_f (trim(fname),H5F_ACC_TRUNC_F,file_id,ierr)
!
! ****** Create the dataspace.
!
      call h5Screate_simple_f (s%ndim,s_dims,dspace_id,ierr)
!
! ****** Create and write the dataset (convert s%f to proper type).
!
      if (s%hdf32) then
        allocate (f4(s_dims(1),s_dims(2),s_dims(3)))
        f4(:,:,:)=s%f(:,:,:)
        call h5Dcreate_f (file_id,'Data',H5T_NATIVE_REAL,dspace_id,dset_id,ierr)
        call h5Dwrite_f (dset_id,H5T_NATIVE_REAL,f4,s_dims,ierr)
        deallocate (f4)
      else
        allocate (f8(s_dims(1),s_dims(2),s_dims(3)))
        f8(:,:,:)=s%f(:,:,:)
        call h5Dcreate_f (file_id,'Data',H5T_NATIVE_DOUBLE,dspace_id,dset_id,ierr)
        call h5Dwrite_f (dset_id,H5T_NATIVE_DOUBLE,f8,s_dims,ierr)
        deallocate (f8)
      endif
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRH5:'
        write (*,*) '### Could not write the dataset.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) '[Error return (from h5Dwrite_f) = ',ierr,']'
        ierr=4
        return
      end if
!
! ****** Check for scales.  If present, add them to the hdf5 dataset.
!
      if (s%scale) then
        do i=1,s%ndim
          if (i.eq.1) then
            dimname='dim1'
          elseif (i.eq.2) then
            dimname='dim2'
          elseif (i.eq.3) then
            dimname='dim3'
          endif
          s_dims_i=s_dims(i)
          call h5Screate_simple_f(1,s_dims_i,dspacedim_id,ierr)
          if (s%hdf32) then
            allocate (f4dim(s_dims_i(1)))
            f4dim(:)=s%scales(i)%f(:)
            call h5Dcreate_f (file_id,dimname,H5T_NATIVE_REAL,dspacedim_id,dim_id,ierr)
            call h5Dwrite_f (dim_id,H5T_NATIVE_REAL,f4dim,s_dims_i,ierr)
            deallocate (f4dim)
          else
            allocate (f8dim(s_dims_i(1)))
            f8dim(:)=s%scales(i)%f(:)
            call h5Dcreate_f (file_id,dimname,H5T_NATIVE_DOUBLE,dspacedim_id,dim_id,ierr)
            call h5Dwrite_f (dim_id,H5T_NATIVE_DOUBLE,f8dim,s_dims_i,ierr)
            deallocate (f8dim)
          endif
          call h5DSset_scale_f (dim_id,ierr,dimname)
          call h5DSattach_scale_f (dset_id,dim_id,i,ierr)
          call h5DSset_label_f(dset_id, i, dimname, ierr)
          call h5Dclose_f (dim_id,ierr)
          call h5Sclose_f (dspacedim_id,ierr)
        end do
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in WRH5:'
          write (*,*) '### Could not write the scales.'
          write (*,*) 'File name: ',trim(fname)
          ierr=5
          return
        endif
      endif
!
! ****** Close the dataset.
!
      call h5Dclose_f (dset_id,ierr)
!
! ****** Close the dataspace.
!
      call h5Sclose_f (dspace_id,ierr)
!
! ****** Close the file.
!
      call h5Fclose_f (file_id,ierr)
!
! ****** Close the hdf5 interface.
!
      call h5close_f (ierr)
!
      end subroutine
!#######################################################################
      subroutine rdtxt_1d (fname,scale,nx,f,x,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read a 1D scientific data set from a text file.
!
!-----------------------------------------------------------------------
!
! ****** This routine calls routine RDTXT to read the file.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      logical :: scale
      integer :: nx
      real(r_typ), dimension(:), pointer :: f
      real(r_typ), dimension(:), pointer :: x
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ierr
!
!-----------------------------------------------------------------------
!
! ****** Declaration for the SDS structure.
!
      type(sds) :: s
!
!-----------------------------------------------------------------------
!
! ****** Read the data set.
!
      call rdtxt (fname,s,ierr)
!
      if (ierr.ne.0) return
!
! ****** Check that this is a 1D data set.
!
      if (s%ndim.ne.1) then
        write (*,*)
        write (*,*) '### ERROR in RDTXT_1D:'
        write (*,*) '### The test file does not contain a 1D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
!
! ****** Set the output arguments.
!
      nx=s%dims(1)
      scale=s%scale
      x=>s%scales(1)%f
!
      allocate (f(nx))
      f=s%f(:,1,1)
      deallocate (s%f)
!
      return
      end
!#######################################################################
      subroutine rdtxt_2d (fname,scale,nx,ny,f,x,y,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read a 2D scientific data set from a text file.
!
!-----------------------------------------------------------------------
!
! ****** This routine calls routine RDTXT to read the file.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      logical :: scale
      integer :: nx,ny
      real(r_typ), dimension(:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: x,y
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ny,ierr
!
!-----------------------------------------------------------------------
!
! ****** Declaration for the SDS structure.
!
      type(sds) :: s
!
!-----------------------------------------------------------------------
!
! ****** Read the data set.
!
      call rdtxt (fname,s,ierr)
!
      if (ierr.ne.0) return
!
! ****** Check that this is a 2D data set.
!
      if (s%ndim.ne.2) then
        write (*,*)
        write (*,*) '### ERROR in RDTXT_2D:'
        write (*,*) '### The text file does not contain a 2D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
!
! ****** Set the output arguments.
!
      nx=s%dims(1)
      ny=s%dims(2)
      scale=s%scale
      x=>s%scales(1)%f
      y=>s%scales(2)%f
!
      allocate (f(nx,ny))
      f=s%f(:,:,1)
      deallocate (s%f)
!
      return
      end
!#######################################################################
      subroutine rdtxt_3d (fname,scale,nx,ny,nz,f,x,y,z,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read a 3D scientific data set from a text file.
!
!-----------------------------------------------------------------------
!
! ****** This routine calls routine RDTXT to read the file.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      logical :: scale
      integer :: nx,ny,nz
      real(r_typ), dimension(:,:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: x,y,z
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ny,nz,ierr
!
!-----------------------------------------------------------------------
!
! ****** Declaration for the SDS structure.
!
      type(sds) :: s
!
!-----------------------------------------------------------------------
!
! ****** Read the data set.
!
      call rdtxt (fname,s,ierr)
!
      if (ierr.ne.0) return
!
! ****** Check that this is a 3D data set.
!
      if (s%ndim.ne.3) then
        write (*,*)
        write (*,*) '### ERROR in RDTXT_3D:'
        write (*,*) '### The text file does not contain a 3D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
!
! ****** Set the output arguments.
!
      nx=s%dims(1)
      ny=s%dims(2)
      nz=s%dims(3)
      scale=s%scale
      x=>s%scales(1)%f
      y=>s%scales(2)%f
      z=>s%scales(3)%f
      f=>s%f
!
      return
      end
!#######################################################################
      subroutine wrtxt_1d (fname,scale,nx,f,x,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Write a 1D scientific data set to a text file.
!
!-----------------------------------------------------------------------
!
! ****** This routine calls routine WRTXT to write the file.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      logical :: scale
      integer :: nx
      real(r_typ), dimension(nx,1,1) :: f
      real(r_typ), dimension(nx) :: x
      integer :: ierr
      intent(in) :: fname,scale,nx,f,x
      intent(out) :: ierr
!
!-----------------------------------------------------------------------
!
! ****** Declaration for the SDS structure.
!
      type(sds) :: s
!
!-----------------------------------------------------------------------
!
! ****** Set the structure components.
!
      s%ndim=1
      s%dims(1)=nx
      s%dims(2)=1
      s%dims(3)=1
      s%scale=scale
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
      else
        nullify (s%scales(1)%f)
      end if
      nullify (s%scales(2)%f)
      nullify (s%scales(3)%f)
      call assign_ptr_3d (f,s%f)
!
! ****** Write the data set.
!
      call wrtxt (fname,s,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRTXT_1D:'
        write (*,*) '### Could not write the 1D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
!
      return
      end
!#######################################################################
      subroutine wrtxt_2d (fname,scale,nx,ny,f,x,y,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Write a 2D scientific data set to a text file.
!
!-----------------------------------------------------------------------
!
! ****** This routine calls routine WRTXT to write the file.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      logical :: scale
      integer :: nx,ny
      real(r_typ), dimension(nx,ny,1) :: f
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      integer :: ierr
      intent(in) :: fname,scale,nx,ny,f,x,y
      intent(out) :: ierr
!
!-----------------------------------------------------------------------
!
! ****** Declaration for the SDS structure.
!
      type(sds) :: s
!
!-----------------------------------------------------------------------
!
! ****** Set the structure components.
!
      s%ndim=2
      s%dims(1)=nx
      s%dims(2)=ny
      s%dims(3)=1
      s%scale=scale
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
        call assign_ptr_1d (y,s%scales(2)%f)
      else
        nullify (s%scales(1)%f)
        nullify (s%scales(2)%f)
      end if
      nullify (s%scales(3)%f)
      call assign_ptr_3d (f,s%f)
!
! ****** Write the data set.
!
      call wrtxt (fname,s,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRTXT_2D:'
        write (*,*) '### Could not write the 2D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
!
      return
      end
!#######################################################################
      subroutine wrtxt_3d (fname,scale,nx,ny,nz,f,x,y,z,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Write a 3D scientific data set to a text file.
!
!-----------------------------------------------------------------------
!
! ****** This routine calls routine WRTXT to write the file.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      logical :: scale
      integer :: nx,ny,nz
      real(r_typ), dimension(nx,ny,nz) :: f
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nz) :: z
      integer :: ierr
      intent(in) :: fname,scale,nx,ny,nz,f,x,y,z
      intent(out) :: ierr
!
!-----------------------------------------------------------------------
!
! ****** Declaration for the SDS structure.
!
      type(sds) :: s
!
!-----------------------------------------------------------------------
!
! ****** Set the structure components.
!
      s%ndim=3
      s%dims(1)=nx
      s%dims(2)=ny
      s%dims(3)=nz
      s%scale=scale
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
        call assign_ptr_1d (y,s%scales(2)%f)
        call assign_ptr_1d (z,s%scales(3)%f)
      else
        nullify (s%scales(1)%f)
        nullify (s%scales(2)%f)
        nullify (s%scales(3)%f)
      end if
      call assign_ptr_3d (f,s%f)
!
! ****** Write the data set.
!
      call wrtxt (fname,s,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRTXT_3D:'
        write (*,*) '### Could not write the 3D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
!
      return
      end
!#######################################################################
      subroutine rdtxt (fname,s,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read a 1D, 2D, or 3D scientific data set from a text file.
!
!-----------------------------------------------------------------------
!
! ****** This routine allocates the required memory and returns
! ****** pointers to the data and scale arrays.
!
!-----------------------------------------------------------------------
!
! ****** Input arguments:
!
!          FNAME   : [character(*)]
!                    Text data file name to read from.
!
! ****** Output arguments:
!
!          S       : [structure of type SDS]
!                    A structure that holds the field, its
!                    dimensions, and the scales, with the
!                    components described below.
!
!          IERR    : [integer]
!                    IERR=0 is returned if the data set was read
!                    successfully.  Otherwise, IERR is set to a
!                    nonzero value.
!
! ****** Components of structure S:
!
!          NDIM    : [integer]
!                    Number of dimensions found in the data set.
!
!          DIMS    : [integer, dimension(3)]
!                    Number of points in the data set dimensions.
!                    For a 1D data set, DIMS(2)=DIMS(3)=1.
!                    For a 2D data set, DIMS(3)=1.
!
!          SCALE   : [logical]
!                    Flag to indicate the presence of scales (axes)
!                    in the data set.  SCALE=.false. means that scales
!                    were not found; SCALE=.true. means that scales
!                    were found.
!
!          HDF32   : [logical]
!                    This flag is is not relevant to text files;
!                    it is used for HDF data.
!                    It is arbitrarily set to HDF32=.false..
!
!          SCALES  : [structure of type RP1D, dimension(3)]
!                    This array holds the pointers to the scales
!                    when SCALE=.true., and is undefined otherwise.
!
!          F       : [real, pointer to a rank-3 array]
!                    This array holds the data set values.
!
! ****** The storage for the arrays pointed to by F, and the
! ****** scales (if present) in structure SCALES, is allocated by
! ****** this routine.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fname
      intent(out) :: s,ierr
!
!-----------------------------------------------------------------------
!
      integer, parameter :: mxdim=3
      integer :: ifscale,i,n
!
!-----------------------------------------------------------------------
!
      ierr=0
!
! ****** Open the file for reading.
!
      call ffopen (1,fname,'r',ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in RDTXT:'
        write (*,*) '### Could not open the text file.'
        write (*,*) 'File name: ',trim(fname)
        ierr=1
        return
      end if
!
! ****** Get the number of dimensions.
!
      call rdint (1,1,s%ndim,ierr)
      if (ierr.ne.0) go to 910
!
      if (s%ndim.lt.1.or.s%ndim.gt.mxdim) then
        write (*,*)
        write (*,*) '### ERROR in RDTXT:'
        write (*,*) '### Invalid number of dimensions in file.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) 'Number of dimensions = ',s%ndim
        write (*,*) 'Maximum number of dimensions = ',mxdim
        ierr=2
        return
      end if
!
! ****** Read the dimensions.
!
      s%dims(:)=1
!
      do i=1,s%ndim
        call rdint (1,1,s%dims(i),ierr)
        if (ierr.ne.0) go to 910
        if (s%dims(i).le.0) go to 920
      enddo
!
! ****** Check if the scales are present.
!
      call rdint (1,1,ifscale,ierr)
      if (ierr.ne.0) go to 910
!
      s%scale=ifscale.ne.0
!
! ****** Allocate memory and read the scales (if present).
!
! ****** If scales are not present,  allocate dummy scales
! ****** (of length 1) so that the pointers to the scales
! ****** are valid.
!
      if (s%scale) then
        do i=1,s%ndim
          allocate (s%scales(i)%f(s%dims(i)))
          call rdfp (1,s%dims(i),s%scales(i)%f,ierr)
          if (ierr.ne.0) go to 910
        enddo
        do i=s%ndim+1,3
          allocate (s%scales(i)%f(1))
        enddo
      else
        allocate (s%scales(1)%f(1))
        allocate (s%scales(2)%f(1))
        allocate (s%scales(3)%f(1))
      end if
!
! ****** Allocate memory for the array.
!
      allocate (s%f(s%dims(1),s%dims(2),s%dims(3)))
!
! ****** Read the data array.
!
      n=product(s%dims(1:s%ndim))
      call rdfp (1,n,s%f,ierr)
      if (ierr.ne.0) go to 910
!
      s%hdf32=.false.
!
      close (1)
!
      return
!
  910 continue
!
      write (*,*)
      write (*,*) '### ERROR in RDTXT:'
      write (*,*) '### Error while reading text data.'
      write (*,*) 'File name: ',trim(fname)
      ierr=3
      return
!
  920 continue
!
      write (*,*)
      write (*,*) '### ERROR in RDTXT:'
      write (*,*) '### Invalid value for dimension.'
      write (*,*) 'Dimension number = ',i
      write (*,*) 'Dimension value = ',s%dims(i)
      ierr=4
!
      return
      end
!#######################################################################
      subroutine rdint (iun,n,i,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read N words of INTEGER data into array I from unit IUN
! ****** using a free format text read.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: iun
      integer :: n
      integer :: i
      integer :: ierr
      intent(in) :: iun,n
      intent(out) :: i,ierr
!
!-----------------------------------------------------------------------
!
      ierr=0
!
      read (iun,*,err=100,end=100) i
!
      return
!
  100 continue
!
! ****** Error in reading the data.
!
      ierr=1
!
      return
      end
!#######################################################################
      subroutine rdfp (iun,n,f,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read N words of REAL data into array F from unit IUN
! ****** using a free format text read.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: iun
      integer :: n
      real(r_typ), dimension(n) :: f
      integer :: ierr
      intent(in) :: iun,n
      intent(out) :: f,ierr
!
!-----------------------------------------------------------------------
!
      ierr=0
!
      read (iun,*,err=100,end=100) f
!
      return
!
  100 continue
!
! ****** Error in reading the data.
!
      ierr=1
!
      return
      end
!#######################################################################
      subroutine wrtxt (fname,s,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Write a 1D, 2D, or 3D scientific data set to a text file.
!
!-----------------------------------------------------------------------
!
! ****** Input arguments:
!
!          FNAME   : [character(*)]
!                    Text data file name to write to.
!
!          S       : [structure of type SDS]
!                    A structure that holds the field, its
!                    dimensions, and the scales, with the
!                    components described below.
!
! ****** Output arguments:
!
!          IERR    : [integer]
!                    IERR=0 is returned if the data set was written
!                    successfully.  Otherwise, IERR is set to a
!                    nonzero value.
!
! ****** Components of structure S:
!
!          NDIM    : [integer]
!                    Number of dimensions in the data set.
!
!          DIMS    : [integer, dimension(3)]
!                    Number of points in the data set dimensions.
!                    Only DIMS(1 .. NDIM) are referenced.
!
!          SCALE   : [logical]
!                    Flag to indicate the presence of scales (axes)
!                    in the data set.  SCALE=.false. means that scales
!                    are not being supplied; SCALE=.true. means that
!                    scales are being supplied.
!
!          HDF32   : [logical]
!                    Flag that indicates the precision of the data.
!                    This flag is used to determine the format for data
!                    written to the text file.  When HDF32=.TRUE., the
!                    data is assumed to originate from a 32-bit HDF data
!                    file, and is written with 7 digits to the text file.
!                    Otherwise, the data is assumed to originate from a
!                    64-bit HDF data file, and is written with 14 digits
!                    to the text file.
!
!          SCALES  : [structure of type RP1D, dimension(3)]
!                    This array holds the pointers to the scales
!                    when SCALE=.true., and is not referenced
!                    otherwise.
!
!          F       : [real, pointer to a rank-3 array]
!                    This array holds the data set values.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fname,s
      intent(out) :: ierr
!
!-----------------------------------------------------------------------
!
! ****** Declarations for temporary variables.
!
      integer :: i
      integer :: n
      character(32) :: fmt
!
!-----------------------------------------------------------------------
!
      ierr=0
!
! ****** Open the file for writing.
!
      call ffopen (1,fname,'rw',ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRTXT:'
        write (*,*) '### Could not open the text file for writing.'
        write (*,*) 'File name: ',trim(fname)
        ierr=1
        return
      end if
!
! ****** Check the number of dimensions.
!
      if (s%ndim.le.0.or.s%ndim.gt.3) then
        write (*,*)
        write (*,*) '### ERROR in WRTXT:'
        write (*,*) '### Could not write the SDS data.'
        write (*,*) 'Invalid number of dimensions.'
        write (*,*) 'NDIM = ',s%ndim
        write (*,*) 'File name: ',trim(fname)
        ierr=1
        return
      end if
!
! ****** Construct the format string for writing floating point
! ****** numbers to the output file.
!
      if (s%hdf32) then
        fmt='(5(1x,1pe13.6))'
      else
        fmt='(3(1x,1pe21.14))'
      end if
!
! ****** Write the number of dimensions.
!
      call wrint (1,1,s%ndim,ierr)
      if (ierr.ne.0) go to 900
!
! ****** Write the dimensions.
!
      do i=1,s%ndim
        call wrint (1,1,s%dims(i),ierr)
        if (ierr.ne.0) go to 900
      enddo
!
! ****** Write the scales.
!
      if (s%scale) then
        call wrint (1,1,1,ierr)
        if (ierr.ne.0) go to 900
        do i=1,s%ndim
          call wrfp (1,s%dims(i),s%scales(i)%f,fmt,ierr)
          if (ierr.ne.0) go to 900
        enddo
      else
        call wrint (1,1,0,ierr)
        if (ierr.ne.0) go to 900
      end if
!
! ****** Write the array.
!
      n=product(s%dims(1:s%ndim))
      call wrfp (1,n,s%f,fmt,ierr)
      if (ierr.ne.0) go to 900
!
      close (1)
!
      return
!
  900 continue
!
      write (*,*)
      write (*,*) '### ERROR in WRTXT:'
      write (*,*) '### Error in writing data to the text file.'
      write (*,*) 'File name: ',trim(fname)
      ierr=2
!
      return
      end
!#######################################################################
      subroutine wrint (iun,n,i,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Write N words of INTEGER data from array I to the file
! ****** connected to unit IUN using a free format text write.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: iun
      integer :: n
      integer :: i
      integer :: ierr
      intent(in) :: iun,n,i
      intent(out) :: ierr
!
!-----------------------------------------------------------------------
!
      ierr=0
!
      write (iun,*,err=100) i
!
      return
!
  100 continue
!
! ****** Error in writing the data.
!
      ierr=1
!
      return
      end
!#######################################################################
      subroutine wrfp (iun,n,f,fmt,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Write N words of REAL data from array F to the file
! ****** connected to unit IUN.
!
! ****** FMT specifies the format string to use.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: iun
      integer :: n
      real(r_typ), dimension(n) :: f
      character(*) :: fmt
      integer :: ierr
      intent(in) :: iun,n,f,fmt
      intent(out) :: ierr
!
!-----------------------------------------------------------------------
!
      ierr=0
!
      write (iun,fmt=fmt,err=100) f
!
      return
!
  100 continue
!
! ****** Error in writing the data.
!
      ierr=1
!
      return
      end
