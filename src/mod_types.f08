module mod_types
  use mod_global_variables
  implicit none

  type ef_type
    integer                  :: index
    character(3)             :: var
    character(str_len)       :: savename
    complex(dp), allocatable :: eigenfunctions(:, :)
  end type ef_type

end module mod_types
