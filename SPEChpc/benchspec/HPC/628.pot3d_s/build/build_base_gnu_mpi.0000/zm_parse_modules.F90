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
      module syntax
!
!-----------------------------------------------------------------------
! ****** Group definitions for parsing command-line arguments.
!-----------------------------------------------------------------------
!
!        GROUP 1: <kw>
!        GROUP 2: <arg>
!        GROUP 3: <kw> <arg>
!        GROUP 4: <kw> <arg> <arg>
!        GROUP 5: GROUP SET (only one must be specified)
!        GROUP 6: GROUP SET (only one or none must be specified)
!
      integer, parameter :: ngroups=6
!
      integer, parameter :: GROUP_K            =1
      integer, parameter :: GROUP_A            =2
      integer, parameter :: GROUP_KA           =3
      integer, parameter :: GROUP_KAA          =4
      integer, parameter :: GROUP_K_ONE_ONLY   =5
      integer, parameter :: GROUP_K_ONE_OR_NONE=6
!
      end module
!#######################################################################
      module string_def
!
!-----------------------------------------------------------------------
! ****** Define a structure to hold a string.
!-----------------------------------------------------------------------
!
      implicit none
!
      type :: string
        character, dimension(:), pointer :: c
      end type
!
      end module
!#######################################################################
      module paragraph_def
!
      use string_def
!
      implicit none
!
!-----------------------------------------------------------------------
! ****** Define a structure for a linked list of lines
! ****** (i.e., a paragraph).
!-----------------------------------------------------------------------
!
      type :: paragraph
        type(string) :: line
        type(paragraph), pointer :: next
      end type
!
!-----------------------------------------------------------------------
! ****** Define a structure to hold a list of paragraphs.
!-----------------------------------------------------------------------
!
      type :: parlist
        type(paragraph), pointer :: par
      end type
!
      end module
!#######################################################################
      module lcase_interface
      interface
        function lcase (s)
        character(*) :: s
        character(len(s)) :: lcase
        end function
      end interface
      end module
!#######################################################################
      module ucase_interface
      interface
        function ucase (s)
        character(*) :: s
        character(len(s)) :: ucase
        end function
      end interface
      end module
!#######################################################################
      module new_par_interface
      interface
        subroutine new_par (par)
        use paragraph_def
        implicit none
        type(paragraph), pointer :: par
        end subroutine
      end interface
      end module
!#######################################################################
      module delete_par_interface
      interface
        subroutine delete_par (par)
        use paragraph_def
        implicit none
        type(paragraph), pointer :: par
        end subroutine
      end interface
      end module
!#######################################################################
      module add_line_interface
      interface
        subroutine add_line (line,par)
        use paragraph_def
        implicit none
        character(*) :: line
        type(paragraph), pointer :: par
        end subroutine
      end interface
      end module
!#######################################################################
      module print_par_interface
      interface
        subroutine print_par (par)
        use paragraph_def
        implicit none
        type(paragraph), pointer :: par
        end subroutine
      end interface
      end module
!#######################################################################
      module get_str_interface
      interface
        function get_str (str)
        use string_def
        implicit none
        type(string) :: str
        character(size(str%c)) :: get_str
        end function
      end interface
      end module
!#######################################################################
      module get_usage_line_interface
      interface
        subroutine get_usage_line (usage)
        use paragraph_def
        implicit none
        type(paragraph), pointer :: usage
        end subroutine
      end interface
      end module
