module mod_types
  use mod_global_variables, only: dp, str_len
  implicit none
  
  private

  type ef_type
    integer                  :: index
    character(3)             :: var
    character(str_len)       :: savename
    complex(dp), allocatable :: eigenfunctions(:, :)
  end type ef_type
  
  public :: ef_type

end module mod_types
