!*******************************************************************
! Fortran interface module to the  Common timer file for use 
! with hpc Weak Scaling              
!
!*******************************************************************
module specmpitime_mod

  use iso_c_binding
  
  interface

  subroutine spectime_start(timer) bind(C)
    use iso_c_binding
    integer(c_int), value :: timer
  end subroutine

  subroutine spectime_stop(timer) bind(C)
    use iso_c_binding
    integer(c_int), value :: timer
  end subroutine

  subroutine spectime_final(pass, units) bind(C)
    use iso_c_binding
    logical,value :: pass 
    real(c_double), value :: units
  end subroutine

  end interface 

end module specmpitime_mod

