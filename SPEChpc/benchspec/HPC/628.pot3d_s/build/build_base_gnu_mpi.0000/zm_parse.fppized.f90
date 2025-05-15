!
!-----------------------------------------------------------------------
!
! ****** Source to build the parsing library.
! ****** These routines are used by Zoran Mikic's tools.
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
!
!-----------------------------------------------------------------------
!
!        07/29/2003, ZM, Version 1.00:
!
!         - Original version of the parsing library.
!           This library was put together to facilitate the
!           development of ZM's tools.
!           It includes routines to parse the command line.
!           The code was cleaned up to use standard FORTRAN90.
!
!        01/17/2005, ZM, Version 1.01:
!
!         - Added the function NARGS_IFIED to return the
!           number of arguments specified on the command line.
!
!        10/28/2005, ZM, Version 1.02:
!
!         - Added the functions LOAD_LIST_OF_REALS and
!           LOAD_LIST_OF_INTS that can be used to parse
!           arguments that contain lists of real or integer values.
!         - Changed the length of temporary arguments to equal
!           512 characters to accomodate long file names.
!
!        10/31/2006, ZM, Version 1.03:
!
!         - Removed the EXTERNAL declarations for GETARG and IARGC
!           since these are now intrinsic routines in the
!           Intel 9.1 compiler.
!
!        03/10/2008, ZM, Version 1.04:
!
!         - Added the LCASE and UCASE functions to convert strings
!           to lowercase and uppercase.
!
!        04/27/2018, ZM, Version 1.05:
!
!         - Added the ability to specify group-sets of keywords.
!           This extends the flexibility of the syntax allowed
!           for keywords.  For example, it allows a syntax in which
!           only one of a set of keywords is allowed to be specified
!           at a time.
!
!        01/08/2019, RC, Version 1.06:
!
!         - Added "append" option to ffopen.
!
!-----------------------------------------------------------------------
!
!#######################################################################
      module parselib_ident
!
      character(*), parameter :: cname='PARSELIB'
      character(*), parameter :: cvers='1.06'
      character(*), parameter :: cdate='01/08/2019'
!
      end module
!#######################################################################
      module parse_args
!
      use string_def
!
      implicit none
!
!-----------------------------------------------------------------------
! ****** Argument descriptor and storage for command-line arguments.
!-----------------------------------------------------------------------
!
! ****** Structure to hold a group-set.
!
      type :: group_set
        integer :: group_type
        integer :: n_members
        integer, dimension(:), allocatable :: members
        logical :: required
        logical :: only_one
      end type
!
! ****** Structure to hold an argument.
!
      type :: arg_descriptor
        integer :: group
        logical :: set
        logical :: required
        type(string) :: keyword
        type(string) :: name
        type(string) :: value
      end type
!
! ****** Maximum number of arguments.
!
      integer, parameter :: mxarg=100
!
! ****** Number of arguments defined.
!
      integer :: nargs
!
! ****** Argument descriptor.
!
      type(arg_descriptor), dimension(mxarg) :: args
!
! ****** Number of arguments specified.
!
      integer :: nargs_spec
!
! ****** Maximum number of group-sets.
!
      integer, parameter :: mxgroup_sets=10
!
! ****** Number of group-sets defined.
!
      integer :: ngroup_sets=0
!
! ****** Group-set descriptor.
!
      type(group_set), dimension(mxgroup_sets) :: group_sets
!
      end module
!#######################################################################
      subroutine ffopen (iun,fname,mode,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Open file FNAME and link it to unit IUN.
!
! ****** If there is an error, this routine returns IERR.ne.0.
!
!-----------------------------------------------------------------------
!
! ****** When MODE='r', the file must exist.
! ****** When MODE='w', the file is created.
! ****** When MODE='rw', the file must exist, but can be overwritten.
! ****** When MODE='a', the file is created if it does not exist,
! ******                otherwise, it is appended.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: iun
      character(*) :: fname
      character(*) :: mode
      integer :: ierr
      logical :: ex
!
!-----------------------------------------------------------------------
!
      ierr=0
!
      if (mode.eq.'r') then
        open (iun,file=fname,status='old',err=900)
      else if (mode.eq.'rw') then
        open (iun,file=fname,status='replace',err=900)
      else if (mode.eq.'w') then
        open (iun,file=fname,form="FORMATTED",status='new',err=900)
      elseif (mode.eq.'a') then
        inquire(file=fname, exist=ex)
        if (ex) then
          open (iun,file=fname,form="FORMATTED",position='append',err=900)
        else
          open (iun,file=fname,form="FORMATTED",status='new',err=900)
        end if
      else
        write (*,*)
        write (*,*) '### ERROR in FFOPEN:'
        write (*,*) '### Invalid MODE requested.'
        write (*,*) 'MODE = ',mode
        write (*,*) 'File name: ',trim(fname)
        ierr=2
        return
      end if
!
      return
!
  900 continue
!
      write (*,*)
      write (*,*) '### ERROR in FFOPEN:'
      write (*,*) '### Error while opening the requested file.'
      write (*,*) 'File name: ',trim(fname)
      write (*,*) 'MODE = ',mode
      ierr=1
!
      return
      end
!#######################################################################
      function lcase (s)
!
!-----------------------------------------------------------------------
!
! ****** Convert the string S into lowercase letters and return it as
! ****** the function result.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*), intent(in) :: s
      character(len(s)) :: lcase
!
!-----------------------------------------------------------------------
!
      integer :: i,ic
!
!-----------------------------------------------------------------------
!
      lcase=' '
!
      do i=1,len_trim(s)
        ic=iachar(s(i:i))
        if (ic.ge.65.and.ic.le.90) then
          ic=ic+32
        end if
        lcase(i:i)=achar(ic)
      end do
!
      return
      end
!#######################################################################
      function ucase (s)
!
!-----------------------------------------------------------------------
!
! ****** Convert the string S into uppercase letters and return it as
! ****** the function result.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*), intent(in) :: s
      character(len(s)) :: ucase
!
!-----------------------------------------------------------------------
!
      integer :: i,ic
!
!-----------------------------------------------------------------------
!
      ucase=' '
!
      do i=1,len_trim(s)
        ic=iachar(s(i:i))
        if (ic.ge.97.and.ic.le.122) then
          ic=ic-32
        end if
        ucase(i:i)=achar(ic)
      end do
!
      return
      end
!#######################################################################
      subroutine parse (errmsg,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Parse the command line.
!
!-----------------------------------------------------------------------
!
! ****** The syntax for the keyword/argument items can be defined
! ****** by using routine DEFARG.
!
! ****** On return, IERR=0 indicates that the command line was
! ****** parsed successfully.
!
! ****** IERR=1 indicates that no arguments were present.  This
! ****** is usually used to print the usage line.
!
! ****** IERR=2 indicates that a syntax error occured.
!
! ****** IERR=3 indicates that one or more required arguments
! ****** was not supplied.
!
! ****** When IERR=2 or IERR=3, an error message is put into
! ****** character string ERRMSG.
!
!-----------------------------------------------------------------------
!
      use syntax
      use parse_args
      use get_str_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: errmsg
      integer :: ierr
!
!-----------------------------------------------------------------------
!
! ****** Command line arguments.
!
      integer :: iargc
!
!-----------------------------------------------------------------------
!
! ****** Carriage return.
!
      character(1), parameter :: cr=achar(10)
!
!-----------------------------------------------------------------------
!
      integer, external :: matchkw
      integer, external :: nextarg
!
!-----------------------------------------------------------------------
!
      character(512) :: arg
      integer :: na,ia,ia0,iarg,ls,i,j,n
!
!-----------------------------------------------------------------------
!
! ****** Initialization.
!
      ierr=0
      nargs_spec=0
      errmsg=' '
!
! ****** Get the number of command line arguments.
!
      na=iargc()
      if (na.eq.0) then
        ierr=1
        return
      end if
!
      ia=1
  200 continue
!
      ia0=ia
!
! ****** Process arguments with syntax: <kw> <arg>
!
      if (na-ia+1.ge.2) then
        call getarg (ia,arg)
        iarg=matchkw(GROUP_KA,trim(arg))
        if (iarg.gt.0) then
          if (.not.args(iarg)%set) then
            ia=ia+1
            call getarg (ia,arg)
            call delete_str (args(iarg)%value)
            call put_str (trim(arg),args(iarg)%value)
            args(iarg)%set=.true.
            ia=ia+1
            nargs_spec=nargs_spec+1
            go to 300
          end if
        end if
      end if
!
! ****** Process arguments with syntax: <kw> <arg> <arg>
!
      if (na-ia+1.ge.3) then
        call getarg (ia,arg)
        iarg=matchkw(GROUP_KAA,trim(arg))
        if (iarg.gt.0) then
          if (.not.args(iarg)%set) then
            ia=ia+1
            call getarg (ia,arg)
            ls=len_trim(arg)
            ls=ls+1
            arg(ls:ls)=' '
            ia=ia+1
            call getarg (ia,arg(ls+1:))
            call delete_str (args(iarg)%value)
            call put_str (trim(arg),args(iarg)%value)
            args(iarg)%set=.true.
            ia=ia+1
            nargs_spec=nargs_spec+1
            go to 300
          end if
        end if
      end if
!
! ****** Process arguments with syntax: <kw>
!
      if (na-ia+1.ge.1) then
        call getarg (ia,arg)
        iarg=matchkw(GROUP_K,trim(arg))
        if (iarg.gt.0) then
          if (.not.args(iarg)%set) then
            call delete_str (args(iarg)%value)
            call put_str (' ',args(iarg)%value)
            args(iarg)%set=.true.
            ia=ia+1
            nargs_spec=nargs_spec+1
            go to 300
          end if
        end if
      end if
!
! ****** Process arguments with syntax: <arg>
!
      if (na-ia+1.ge.1) then
        iarg=nextarg(GROUP_A)
        if (iarg.gt.0) then
          call getarg (ia,arg)
          call delete_str (args(iarg)%value)
          call put_str (trim(arg),args(iarg)%value)
          args(iarg)%set=.true.
          ia=ia+1
          nargs_spec=nargs_spec+1
          go to 300
        end if
      end if
!
  300 continue
!
! ****** Check that an argument was found.
!
      if (ia.eq.ia0) then
        ierr=2
        errmsg='### Syntax error.'
        return
      end if
!
! ****** Keep processing arguments until done.
!
      if (na-ia+1.gt.0) go to 200
!
! ****** Check that the required arguments were supplied.
!
      do i=1,nargs
        if (args(i)%required.and..not.args(i)%set) then
          ierr=3
          errmsg='### A required argument was not supplied.'
          return
        end if
      enddo
!
! ****** Check that the group-set keywords were specified
! ****** properly.
!
! ****** Check the "only one" and "required" group-set keywords.
!
      do i=1,ngroup_sets
        n=0
        do j=1,group_sets(i)%n_members
          iarg=group_sets(i)%members(j)
          if (args(iarg)%set) n=n+1
        enddo
        if (n.gt.1.and.group_sets(i)%only_one) then
          ierr=2
          errmsg='### Syntax error.'
          errmsg=trim(errmsg)//cr//cr
          errmsg=trim(errmsg)//' ### You can only specify'// &
                ' one of the following keywords:'
          do j=1,group_sets(i)%n_members
            iarg=group_sets(i)%members(j)
            errmsg=trim(errmsg)//cr//' '// &
                  get_str(args(iarg)%keyword)
          enddo
          return
        else if (n.eq.0.and.group_sets(i)%required) then
          ierr=2
          errmsg='### Syntax error.'
          errmsg=trim(errmsg)//cr//cr
          errmsg=trim(errmsg)//' ### You must specify'// &
                ' one of the following keywords:'
          do j=1,group_sets(i)%n_members
            iarg=group_sets(i)%members(j)
            errmsg=trim(errmsg)//cr//' '//  &
                  get_str(args(iarg)%keyword)
          enddo
          return
        end if
      enddo
!
      return
      end
!#######################################################################
      subroutine get_usage_line (usage)
!
!-----------------------------------------------------------------------
!
! ****** Construct the usage line in paragraph USAGE.
!
! ****** Use routine PRINT_PAR to write the usage line.
!
!-----------------------------------------------------------------------
!
      use syntax
      use parse_args
      use paragraph_def
      use new_par_interface
      use add_line_interface
      use get_str_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(paragraph), pointer :: usage
!
!-----------------------------------------------------------------------
!
! ****** Right margin for printing the usage line.
!
      integer, parameter :: rmargin=78
!
!-----------------------------------------------------------------------
!
      character(512) :: line
      integer :: iarg,n0
      type(paragraph), pointer :: current_par
!
      logical, dimension(:), allocatable :: done
      integer :: i,j,ix
      logical :: belongs_to_gs
      character(512) :: str
!
!-----------------------------------------------------------------------
!
! ****** Construct the usage line in USAGE.
!
      call new_par (usage)
      current_par=>usage
!
! ****** Start with the command name (as invoked).
!
      call getarg (0,line)
!
! ****** Allocate storage.
!
      allocate (done(nargs))
      done=.false.
!
      iarg=1
!
! ****** Add the arguments one by one.
!
      do while (iarg.le.nargs)

        if (done(iarg)) then
          iarg=iarg+1
          cycle
        end if
!
! ****** Add the syntax for the next argument to LINE.
!
        n0=len_trim(line)
!
! ****** Check if this argument belongs to a "only one" group-set.
! ****** If so, collect all members of the group-set.  Otherwise, just
! ****** process this argument.
!
        belongs_to_gs=.false.
        do i=1,ngroup_sets
          if (group_sets(i)%only_one) then
            do j=1,group_sets(i)%n_members
              if (group_sets(i)%members(j).eq.iarg) then
                ix=i
                belongs_to_gs=.true.
                exit
              end if
            enddo
          end if
          if (belongs_to_gs) exit
        enddo
!
        if (belongs_to_gs) then
!
! ****** This argument belongs to a "only one" group-set:
! ****** process all its members.
!
          str=''
          do j=1,group_sets(ix)%n_members
            i=group_sets(ix)%members(j)
            done(i)=.true.
            str=trim(str)//get_str(args(i)%keyword)
            if (j.lt.group_sets(ix)%n_members) then
              str=trim(str)//'|'
            end if
          enddo
          if (group_sets(ix)%required) then
            line=trim(line)//' '//trim(str)
          else
            line=trim(line)//' ['//trim(str)//']'
          end if
!
        else
!
! ****** This argument does not belong to a "only one" group-set:
! ****** process just this single argument.
!
          if (args(iarg)%required) then
            line=trim(line)//' '//get_str(args(iarg)%keyword)
          else
            line=trim(line)//' ['//get_str(args(iarg)%keyword)
          end if
          line=trim(line)//' '//get_str(args(iarg)%name)
          if (.not.args(iarg)%required) then
            line=trim(line)//']'
          end if
          done(iarg)=.true.
!
        end if
!
! ****** Check if the addition of the argument causes the line
! ****** to wrap; if it does, break the line prior to the
! ****** argument text.
!
        if (len_trim(line).gt.rmargin) then
          call add_line (line(1:n0),current_par)
          line=' '//line(n0+1:)
        end if
!
! ****** If the line is still too long, force a break at RMARGIN
! ****** until the line is shorter than RMARGIN.
!
        do while (len_trim(line).gt.rmargin)
          call add_line (line(1:rmargin),current_par)
          line='  '//line(rmargin+1:)
        enddo
!
! ****** Process the next argument.
!
        iarg=iarg+1
!
      enddo
!
! ****** Add the last line to the paragraph.
!
      if (line.ne.' ') call add_line (trim(line),current_par)
!
! ****** Deallocate storage.
!
      deallocate (done)
!
      return
      end
!#######################################################################
      subroutine defarg (group,keyword,default,name)
!
!-----------------------------------------------------------------------
!
! ****** Define the syntax for a command line argument item.
!
!-----------------------------------------------------------------------
!
! ****** GROUP is the syntax group index;
! ****** KEYWORD is the keyword;
! ****** DEFAULT is the default value of the argument;
! ****** NAME is the name of the argument (for use in
! ****** constructing the usage line).
!
!-----------------------------------------------------------------------
!
      use syntax
      use parse_args
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: group
      character(*) :: keyword
      character(*) :: default
      character(*) :: name
!
!-----------------------------------------------------------------------
!
      integer :: i,n,ix
      character(512), dimension(:), allocatable :: names
!
!-----------------------------------------------------------------------
!
      integer, external :: number_of_names
      integer, external :: match_keyword
!
!-----------------------------------------------------------------------
!
! ****** Check that the group index is valid.
!
      if (group.lt.0.or.group.gt.ngroups) then
        write (*,*)
        write (*,*) '### ERROR in DEFARG:'
        write (*,*) '### An invalid group index was specified.'
        write (*,*) 'Group index = ',group
        write (*,*) 'Keyword: ',trim(keyword)
        if (name.ne.' ') write (*,*) 'Name: ',trim(name)
        write (*,*)
        write (*,*) '### This indicates a programming error'// &
                   ' in the syntax definition and use.'
        call exit (1)
      end if
!
! ****** Check for a null keyword.
!
      if (keyword.eq.' ') then
        write (*,*)
        write (*,*) '### ERROR in DEFARG:'
        write (*,*) '### The keyword is null.'
        write (*,*) 'Group index = ',group
        if (name.ne.' ') write (*,*) 'Name: ',trim(name)
        write (*,*)
        write (*,*) '### This indicates a programming error'// &
                    ' in the syntax definition and use.'
        call exit (1)
      end if
!
! ****** Check if this defines a group-set.
!
      if (group.eq.GROUP_K_ONE_ONLY.or. &
          group.eq.GROUP_K_ONE_OR_NONE) then
!
! ****** Increment the group-set counter.
!
        if (ngroup_sets.ge.mxgroup_sets) then
          write (*,*)
          write (*,*) '### ERROR in DEFARG:'
          write (*,*) '### Exceeded the number of allowed group-sets.'
          write (*,*) 'Maximum number of group-sets = ',mxgroup_sets
          write (*,*) 'Group index = ',group
          write (*,*) 'Keyword: ',trim(keyword)
          write (*,*)
          write (*,*) '### This indicates a programming error'// &
                      ' in the syntax definition and use.'
          call exit (1)
        end if
!
        ngroup_sets=ngroup_sets+1
!
! ****** Set the group type.
!
        group_sets(ngroup_sets)%group_type=group
!
! ****** Set the appropriate group-set flags.
!
        if (group.eq.GROUP_K_ONE_ONLY) then
          group_sets(ngroup_sets)%required=.true.
          group_sets(ngroup_sets)%only_one=.true.
        else if (group.eq.GROUP_K_ONE_OR_NONE) then
          group_sets(ngroup_sets)%required=.false.
          group_sets(ngroup_sets)%only_one=.true.
        else
          group_sets(ngroup_sets)%required=.false.
          group_sets(ngroup_sets)%only_one=.false.
        end if
!
! ****** Mark the members in the group-set.
!
        n=number_of_names(keyword)
        allocate (names(n))
        allocate (group_sets(ngroup_sets)%members(n))
        group_sets(ngroup_sets)%n_members=n
        call load_list_of_names (keyword,n,names)
        do i=1,n
          ix=match_keyword(names(i))
          if (ix.eq.0) then
            write (*,*)
            write (*,*) '### ERROR in DEFARG:'
            write (*,*) '### Error while defining a group set.'
            write (*,*) 'Could not match a group-set keyword to'// &
                        ' an existing keyword:'
            write (*,*) 'Group-set keyword: ',trim(keyword)
            write (*,*) 'Keyword not matched: ',trim(names(i))
            write (*,*)
            write (*,*) '### This indicates a programming error'// &
                        ' in the syntax definition and use.'
            call exit (1)
          end if
          group_sets(ngroup_sets)%members(i)=ix
        enddo
        deallocate (names)
!
        return
!
      end if
!
! ****** Increment the argument counter.
!
      if (nargs.ge.mxarg) then
        write (*,*)
        write (*,*) '### ERROR in DEFARG:'
        write (*,*) '### Exceeded the number of allowed arguments.'
        write (*,*) 'Maximum number of arguments = ',mxarg
        write (*,*) 'Group index = ',group
        write (*,*) 'Keyword: ',trim(keyword)
        if (name.ne.' ') write (*,*) 'Name: ',trim(name)
        write (*,*)
        write (*,*) '### This indicates a programming error'// &
                    ' in the syntax definition and use.'
        call exit (1)
      end if
!
      nargs=nargs+1
!
! ****** Store the group index and keyword.
!
! ****** For group GROUP_A (single arguments), the name of the
! ****** argument is passed as the "keyword".
!
      args(nargs)%group=group
      call put_str (trim(keyword),args(nargs)%keyword)
!
! ****** Initialize the flag that indicates whether an argument
! ****** has been set.
!
      args(nargs)%set=.false.
!
! ****** If a default argument was supplied, the argument
! ****** does not have to be set.  Use DEFAULT=' ' to
! ****** indicate that an argument is required.
!
! ****** If a default argument has been supplied, store it in
! ****** ARGS(nargs)%VALUE.  If there is no default,
! ****** set ARGS(nargs)%VALUE to an empty string.
!
! ****** Since group GROUP_K doesn't have an argument,
! ****** DEFAULT is ignored for this group.
!
      if (group.eq.GROUP_K) then
        args(nargs)%required=.false.
        call put_str (' ',args(nargs)%value)
      else
        if (default.eq.' ') then
          args(nargs)%required=.true.
          call put_str (' ',args(nargs)%value)
        else
          args(nargs)%required=.false.
          call put_str (trim(default),args(nargs)%value)
        end if
      end if
!
! ****** Store the argument name.  For groups GROUP_K (keywords)
! ****** and GROUP_A (single arguments), there is no argument name,
! ****** so NAME is ignored.
!
      if (group.eq.GROUP_K.or.group.eq.GROUP_A) then
        call put_str (' ',args(nargs)%name)
      else
        call put_str (trim(name),args(nargs)%name)
      end if
!
      return
      end
!#######################################################################
      subroutine fetcharg (keyword,set,arg)
!
!-----------------------------------------------------------------------
!
! ****** Fetch the value of the argument corresponding to
! ****** keyword KEYWORD.
!
!-----------------------------------------------------------------------
!
! ****** If KEYWORD is a keyword-type argument (GROUP_K), return
! ****** its setting through variable SET.  The variable ARG should
! ****** be ignored for this type of keyword.
!
! ****** For keywords with arguments (GROUP_A, GROUP_KA, and
! ****** GROUP_KAA), return the value of the arguments in ARG,
! ****** and return SET=.true. if they were set via the command line;
! ****** otherwise, return SET=.false..
!
!-----------------------------------------------------------------------
!
      use parse_args
      use get_str_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: keyword
      logical :: set
      character(*) :: arg
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
      do i=nargs,1,-1
        if (keyword.eq.get_str(args(i)%keyword)) go to 100
      enddo
!
      write (*,*)
      write (*,*) '### ERROR in FETCHARG:'
      write (*,*) '### The requested keyword could not be matched.'
      write (*,*) 'Keyword = ',trim(keyword)
      write (*,*)
      write (*,*) '### This indicates a programming error'// &
                  ' in the syntax definition and use.'
      call exit (1)
!
  100 continue
!
      set=args(i)%set
      arg=get_str(args(i)%value)
!
      return
      end
!#######################################################################
      function matchkw (group,keyword)
!
!-----------------------------------------------------------------------
!
! ****** Match keyword KEYWORD against the list of keywords in
! ****** group GROUP.
!
! ****** If found, set the function value to the corresponding
! ****** argument number.  Otherwise, return MATCHKW=0.
!
!-----------------------------------------------------------------------
!
      use parse_args
      use get_str_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: group
      character(*) :: keyword
      integer :: matchkw
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
      matchkw=0
!
      do i=nargs,1,-1
        if (group.eq.args(i)%group) then
          if (keyword.eq.get_str(args(i)%keyword)) then
            matchkw=i
            return
          end if
        end if
      enddo
!
      return
      end
!#######################################################################
      function nextarg (group)
!
!-----------------------------------------------------------------------
!
! ****** Find the position of the next argument in group GROUP
! ****** that has not been set.
!
!-----------------------------------------------------------------------
!
! ****** If an empty slot is found, set the function value
! ****** to the corresponding argument number.
!
! ****** Otherwise, return NXTARG=0.
!
!-----------------------------------------------------------------------
!
      use parse_args
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: group
      integer :: nextarg
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
      nextarg=0
!
      do i=1,nargs
        if (group.eq.args(i)%group) then
          if (.not.args(i)%set) then
            nextarg=i
            return
          end if
        end if
      enddo
!
      return
      end
!#######################################################################
      subroutine nargs_specified (n)
!
!-----------------------------------------------------------------------
!
! ****** Return the number of arguments specified on the command
! ****** line.
!
!-----------------------------------------------------------------------
!
      use parse_args
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: n
!
!-----------------------------------------------------------------------
!
      n=nargs_spec
!
      return
      end
!#######################################################################
      function match_keyword (keyword)
!
!-----------------------------------------------------------------------
!
! ****** Match KEYWORD to the current list of defined keywords.
!
! ****** A successful match returns the index of the first matched
! ****** entry in TABLE; an unsuccessful match returns 0.
!
!-----------------------------------------------------------------------
!
      use parse_args
      use get_str_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*), intent(in) :: keyword
      integer :: match_keyword
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
      match_keyword=0
!
      do i=1,nargs
        if (keyword.eq.get_str(args(i)%keyword)) then
          match_keyword=i
          return
        end if
      end do
!
      return
      end
!#######################################################################
      subroutine new_par (par)
!
!-----------------------------------------------------------------------
!
! ****** Initialize paragraph PAR.
!
!-----------------------------------------------------------------------
!
      use paragraph_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(paragraph), pointer :: par
!
!-----------------------------------------------------------------------
!
      allocate (par)
      nullify (par%line%c)
      nullify (par%next)
!
      return
      end
!#######################################################################
      subroutine delete_par (par)
!
!-----------------------------------------------------------------------
!
! ****** Delete paragraph PAR and deallocate its storage and that
! ****** of its linked lists.
!
!-----------------------------------------------------------------------
!
      use paragraph_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(paragraph), pointer :: par
!
!-----------------------------------------------------------------------
!
      type(paragraph), pointer :: current_par,previous_par
!
!-----------------------------------------------------------------------
!
      current_par=>par
!
      do
!
! ****** Deallocate the line buffer.
!
        call delete_str (current_par%line)
!
! ****** Set the pointer to the next line (if it has been defined).
!
        if (.not.associated(current_par%next)) exit
        previous_par=>current_par
        current_par=>current_par%next
        deallocate (previous_par)
!
      enddo
!
      deallocate (current_par)
!
      return
      end
!#######################################################################
      subroutine add_line (line,par)
!
!-----------------------------------------------------------------------
!
! ****** Add LINE to paragraph PAR.
!
! ****** On exit from this routine, PAR points to a new line,
! ****** and can be used to store the next line of text.
!
!-----------------------------------------------------------------------
!
      use paragraph_def
      use new_par_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: line
      type(paragraph), pointer :: par
!
!-----------------------------------------------------------------------
!
! ****** Store LINE into the string buffer for the current line.
!
      call put_str (line,par%line)
!
! ****** Allocate a pointer to the next line.
!
      call new_par (par%next)
!
! ****** Set PAR to point to the next line.
!
      par=>par%next
!
      return
      end
!#######################################################################
      subroutine print_par (par)
!
!-----------------------------------------------------------------------
!
! ****** Print all lines of paragraph PAR to STDOUT.
!
!-----------------------------------------------------------------------
!
      use paragraph_def
      use get_str_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(paragraph), pointer :: par
!
!-----------------------------------------------------------------------
!
      type(paragraph), pointer :: current_par
!
!-----------------------------------------------------------------------
!
      current_par=>par
!
      do
!
! ****** Print the line if it has been defined.
!
        if (associated(current_par%line%c)) then
          write (*,*) trim(get_str(current_par%line))
        end if
!
! ****** Set the pointer to the next line (if it has been defined).
!
        if (.not.associated(current_par%next)) exit
        current_par=>current_par%next
!
      enddo
!
      return
      end
!#######################################################################
      subroutine put_str (cval,str)
!
!-----------------------------------------------------------------------
!
! ****** Store character variable CVAL into string STR.
! ****** This routine allocates storage for the string.
!
!-----------------------------------------------------------------------
!
      use string_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: cval
      type(string) :: str
!
!-----------------------------------------------------------------------
!
      integer :: l,i
!
!-----------------------------------------------------------------------
!
      l=len(cval)
!
      allocate (str%c(l))
!
      do i=1,l
        str%c(i)=cval(i:i)
      enddo
!
      return
      end
!#######################################################################
      function get_str (str)
!
!-----------------------------------------------------------------------
!
! ****** Return the value of string STR as the function value
! ****** (as an assumed-length character variable).
!
!-----------------------------------------------------------------------
!
      use string_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(string) :: str
      character(size(str%c)) :: get_str
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
      do i=1,size(str%c)
        get_str(i:i)=str%c(i)
      enddo
!
      return
      end
!#######################################################################
      subroutine delete_str (str)
!
!-----------------------------------------------------------------------
!
! ****** Delete the storage for string STR.
!
!-----------------------------------------------------------------------
!
      use string_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(string) :: str
!
!-----------------------------------------------------------------------
!
      if (associated(str%c)) then
        deallocate (str%c)
      end if
      nullify (str%c)
!
      return
      end
!#######################################################################
      function intval (avalue,name)
!
!-----------------------------------------------------------------------
!
! ****** Get the value of the integer in character variable AVALUE.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: avalue
      character(*) :: name
      integer :: intval
!
!-----------------------------------------------------------------------
!
      logical, external :: ifint
!
!-----------------------------------------------------------------------
!
      integer :: ivalue
!
!-----------------------------------------------------------------------
!
      if (.not.ifint(trim(avalue),ivalue)) then
        write (*,*)
        write (*,*) '### ERROR in INTVAL:'
        write (*,*) '### Could not interpret an integer '// &
                    'while setting: ',trim(name)
        write (*,*) 'Invalid format: ',trim(avalue)
        call exit (1)
      end if
!
      intval=ivalue
!
      return
      end
!#######################################################################
      function fpval (avalue,name)
!
!-----------------------------------------------------------------------
!
! ****** Get the value of the floating point number in character
! ****** variable AVALUE.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      character(*) :: avalue
      character(*) :: name
      real(r_typ) :: fpval
!
!-----------------------------------------------------------------------
!
      logical, external :: iffp
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: value
!
!-----------------------------------------------------------------------
!
      if (.not.iffp(trim(avalue),value)) then
        write (*,*)
        write (*,*) '### ERROR in FPVAL:'
        write (*,*) '### Could not interpret a floating point '// &
                    'number while setting: ',trim(name)
        write (*,*) 'Invalid format: ',trim(avalue)
        call exit (1)
      end if
!
      fpval=value
!
      return
      end
!#######################################################################
      function iffp (alpha,value)
!
!-----------------------------------------------------------------------
!
! ****** Determine if ALPHA represents a floating point number;
! ****** if so, return its value in VALUE.
!
!-----------------------------------------------------------------------
!
! ****** Set IFFP=.TRUE. if ALPHA contains an alphanumeric
! ****** string with the following format:
!
!       ALPHA = '[A][B...B][.][B...B][e[A]B[B...B]]',
!
! ****** where A represents a + or - sign, and B represents a digit
! ****** between 0 and 9, inclusive.
! ****** The exponent may be denoted by a lower or upper case e.
! ****** The mantissa must have at least one digit, and the
! ****** the exponent, if present, must have between 1 and 3 digits.
! ****** Otherwise, set IFFP=.FALSE.
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
      character(*) :: alpha
      real(r_typ) :: value
      logical :: iffp
!
!-----------------------------------------------------------------------
!
      integer :: nmant,nexp,k,i,ke
      logical :: ifpoint,ifexp
      character(7) :: fmt
!
!-----------------------------------------------------------------------
!
      iffp=.false.
      ifpoint=.false.
      ifexp=.false.
      nmant=0
      nexp=0
!
      do k=1,len_trim(alpha)
        i=iachar(alpha(k:k))
!
! ****** Check for a sign in the first position.
!
        if (k.eq.1.and.(i.eq.43.or.i.eq.45)) cycle
!
! ****** Check for a digit.
!
        if (i.ge.48.and.i.le.57) then
!
! ****** Count digits in mantissa and exponent.
!
        if (ifexp) then
          nexp=nexp+1
          else
            nmant=nmant+1
          end if
          cycle
!
        end if
!
! ****** Check for a decimal point.
!
        if (.not.ifpoint.and.i.eq.46) then
!
! ****** Check that we are in the mantissa.
!
          if (.not.ifexp) then
            ifpoint=.true.
            cycle
          end if
!
        end if
!
! ****** Check for an exponent.
!
        if (.not.ifexp.and.(i.eq.101.or.i.eq.69)) then
          ifexp=.true.
          ke=k
          cycle
        end if
!
! ****** Check for an exponent sign.
!
        if (ifexp.and.k.eq.ke+1.and.(i.eq.43.or.i.eq.45)) cycle
!
! ****** Failed check: fall through here.
!
        iffp=.false.
!
        return
!
      enddo
!
! ****** Final check of validity: check number of digits in
! ****** the mantissa and exponent.
!
      if (nmant.ge.1) iffp=.true.
      if (ifexp.and.(nexp.lt.1.or.nexp.gt.3)) iffp=.false.
!
! ****** Obtain its numeric value.
!
      fmt='(f  .0)'
      write (fmt(3:4),'(i2.2)') len_trim(alpha)
!
      if (iffp) read (alpha,fmt) value
!
      return
      end
!#######################################################################
      function ifint (alpha,ivalue)
!
!-----------------------------------------------------------------------
!
! ****** If ALPHA represents an integer, return IFINT=.true., and
! ****** put its value into IVALUE.
!
! ****** Otherwise, return IFINT=.false..
!
!-----------------------------------------------------------------------
!
! ****** A valid integer has the format:
!
!          ALPHA = '[A]B[B...B]',
!
! ****** where A represents a + or - sign, and B represents a digit
! ****** between 0 and 9, inclusive.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: alpha
      integer :: ivalue
      logical :: ifint
!
!-----------------------------------------------------------------------
!
      integer :: k,i
      character(5) :: fmt
!
!-----------------------------------------------------------------------
!
      ifint=.false.
!
      do k=1,len_trim(alpha)
!
        i=iachar(alpha(k:k))
!
! ****** Check for a sign in the first position.
!
        if (k.eq.1.and.(i.eq.43.or.i.eq.45)) cycle
!
! ****** Check for a digit.
!
        if (i.ge.48.and.i.le.57) then
          ifint=.true.
          cycle
        end if
!
! ****** Failed check: fall through here.
!
        ifint=.false.
!
        return
!
      enddo
!
! ****** Obtain its numeric value.
!
      fmt='(i  )'
      write (fmt(3:4),'(i2.2)') len_trim(alpha)
!
      if (ifint) read (alpha,fmt) ivalue
!
      return
      end
!#######################################################################
      subroutine load_list_of_reals (s,label,n,f)
!
!-----------------------------------------------------------------------
!
! ****** Read N real values from character string S into
! ****** array F(N). The values in S may be either space or
! ****** comma separated.
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
      character(*) :: s
      character(*) :: label
      integer :: n
      real(r_typ), dimension(n) :: f
      intent(in) :: s,label,n
      intent(out) :: f
!
!-----------------------------------------------------------------------
!
      integer :: i,i0,i1
      character :: delimiter
      character(512) :: list
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: fpval
!
!-----------------------------------------------------------------------
!
! ****** Make a local copy of the string (removing leading spaces).
!
      list=adjustl(s)
!
! ****** If any commas are present, use a comma as the delimiter.
! ****** Otherwise, one or more spaces is used as a delimiter.
! ****** In this case, compress multiple spaces into a single space.
!
      if (index(list,',').ne.0) then
        delimiter=','
      else
        delimiter=' '
        call delete_repeated_char (list,' ')
      end if
!
! ****** Read the list of N numbers sequentially into F.
!
      i0=1
      do i=1,n-1
        i1=scan(list(i0:),delimiter)+i0-2
        f(i)=fpval(adjustl(list(i0:i1)),label)
        i0=i1+2
      enddo
      f(n)=fpval(adjustl(list(i0:)),label)
!
      return
      end
!#######################################################################
      subroutine load_list_of_ints (s,label,n,j)
!
!-----------------------------------------------------------------------
!
! ****** Read N integer values from character string S into
! ****** array J(N).  The values in S may be either space or
! ****** comma separated.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: s
      character(*) :: label
      integer :: n
      integer, dimension(n) :: j
      intent(in) :: s,label,n
      intent(out) :: j
!
!-----------------------------------------------------------------------
!
      integer :: i,i0,i1
      character :: delimiter
      character(512) :: list
!
!-----------------------------------------------------------------------
!
      integer, external :: intval
!
!-----------------------------------------------------------------------
!
! ****** Make a local copy of the string (removing leading spaces).
!
      list=adjustl(s)
!
! ****** If any commas are present, use a comma as the delimiter.
! ****** Otherwise, one or more spaces is used as a delimiter.
! ****** In this case, compress multiple spaces into a single space.
!
      if (index(list,',').ne.0) then
        delimiter=','
      else
        delimiter=' '
        call delete_repeated_char (list,' ')
      end if
!
! ****** Read the list of N numbers sequentially into J.
!
      i0=1
      do i=1,n-1
        i1=scan(list(i0:),delimiter)+i0-2
        j(i)=intval(adjustl(list(i0:i1)),label)
        i0=i1+2
      enddo
      j(n)=intval(adjustl(list(i0:)),label)
!
      return
      end
!#######################################################################
      function number_of_names (s)
!
!-----------------------------------------------------------------------
!
! ****** Return the number of names in string S.  The values
! ****** can be separated by either spaces or commas.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: s
      integer :: number_of_names
      intent(in) :: s
!
!-----------------------------------------------------------------------
!
      integer :: l,n,i0,i1
      character :: delimiter
      character(len(s)) :: list
!
!-----------------------------------------------------------------------
!
! ****** Make a local copy of the string (removing leading spaces).
!
      list=adjustl(s)
!
! ****** If any commas are present, use a comma as the delimiter.
! ****** Otherwise, one or more spaces is used as a delimiter.
! ****** In this case, compress multiple spaces into a single space.
!
      if (index(list,',').ne.0) then
        delimiter=','
      else
        delimiter=' '
        call delete_repeated_char (list,' ')
      end if
!
! ****** Find the number of names N.
!
      l=len_trim(list)
!
      n=0
!
      i0=1
      do
        i1=scan(list(i0:l),delimiter)
        if (i1.eq.0) then
          if (.not.(delimiter.eq.' '.and.l.lt.i0)) then
            n=n+1
          end if
          exit
        end if
        i1=i1+i0-2
        i0=i1+2
        n=n+1
      enddo
!
      number_of_names=n
!
      return
      end
!#######################################################################
      subroutine load_list_of_names (s,n,names)
!
!-----------------------------------------------------------------------
!
! ****** Read up to N names from character string S into
! ****** array NAMES(N).  The values can be separated by either
! ****** spaces or commas.
!
!-----------------------------------------------------------------------
!
! ****** This routine is designed to be used in conjunction with
! ****** routine NUMBER_OF_NAMES.  A typical use would be:
!
!          character(32) :: s
!          character(32), dimension(:), allocatable :: names
!          integer :: n
!
!          n=number_of_names(s)
!          allocate (names(n))
!          call load_list_of_names (s,n,names)
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: s
      integer :: n
      character(*), dimension(n) :: names
      intent(in) :: s,n
      intent(out) :: names
!
!-----------------------------------------------------------------------
!
      integer :: l,i,i0,i1
      character :: delimiter
      character(len(s)) :: list
!
!-----------------------------------------------------------------------
!
! ****** Make a local copy of the string (removing leading spaces).
!
      list=adjustl(s)
!
! ****** If any commas are present, use a comma as the delimiter.
! ****** Otherwise, one or more spaces is used as a delimiter.
! ****** In this case, compress multiple spaces into a single space.
!
      if (index(list,',').ne.0) then
        delimiter=','
      else
        delimiter=' '
        call delete_repeated_char (list,' ')
      end if
!
! ****** Read the names.
!
      l=len_trim(list)
!
      i=0
!
      i0=1
      do
        i1=scan(list(i0:l),delimiter)
        if (i1.eq.0) then
          if (.not.(delimiter.eq.' '.and.l.lt.i0)) then
            i=i+1
            if (i.gt.n) exit
            names(i)=adjustl(list(i0:l))
          end if
          exit
        end if
        i=i+1
        if (i.gt.n) exit
        i1=i1+i0-2
        names(i)=adjustl(list(i0:i1))
        i0=i1+2
      enddo
!
      return
      end
!#######################################################################
      subroutine delete_repeated_char (s,c)
!
!-----------------------------------------------------------------------
!
! ****** Transform repeated adjoining occurrences of character C
! ****** in string S into single occurrences of C.
!
! ****** The string S is overwritten by the modified string.
!
! ****** Trailing blanks in S are ignored.
!
!-----------------------------------------------------------------------
!
! ****** For example, suppose this routine is called with C='d' and
! ****** S='abcdddeefdhdd'.  On return, S will have the value
! ****** 'abcdeefdhd'.
!
!-----------------------------------------------------------------------
!
! ****** This routine uses the FORTRAN90 intrinsic SCAN.
!
!-----------------------------------------------------------------------
!
      character(*) :: s
      character :: c
      intent(in) :: c
      intent(inout) :: s
!
!-----------------------------------------------------------------------
!
      integer :: i,i0
!
!-----------------------------------------------------------------------
!
      i0=1
      do
        i=scan(trim(s(i0:)),c)
        if (i.eq.0) exit
        i0=i0+i
        do
          if (s(i0:i0).ne.c) exit
          s(i0:)=s(i0+1:)
        enddo
      enddo
!
      return
      end
