module mod_types
  use mod_global_variables
  implicit none

  type eigenf_type
    integer                  :: index
    integer                  :: write_out
    character(3)             :: var
    character(str_len)       :: savename
    complex(dp), allocatable :: eigenfunctions(:, :)
  end type eigenf_type

end module mod_types
