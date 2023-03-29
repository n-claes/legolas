module mod_gravity
  use mod_global_variables, only: dp
  use mod_function_utils, only: zero_func
  implicit none

  private

  type, public :: gravity_t
    procedure(real(dp)), pointer, nopass :: g0
  contains
    procedure :: delete
  end type gravity_t

  public :: new_gravity

contains

  function new_gravity() result(gravity)
    type(gravity_t) :: gravity
    gravity%g0 => zero_func
  end function new_gravity


  pure subroutine delete(this)
    class(gravity_t), intent(inout) :: this
    nullify(this%g0)
  end subroutine delete

end module mod_gravity
