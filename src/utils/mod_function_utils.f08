module mod_function_utils
  use mod_global_variables, only: dp
  implicit none

  private

  public :: zero_func
  public :: from_function

contains

  real(dp) function zero_func()
    zero_func = 0.0_dp
  end function zero_func


  function from_function(func, values) result(array)
    procedure(real(dp)), pointer :: func
    real(dp), intent(in) :: values(:)
    real(dp) :: array(size(values))
    integer :: i

    do i = 1, size(values)
      array(i) = func(values(i))
    end do
  end function from_function

end module mod_function_utils
