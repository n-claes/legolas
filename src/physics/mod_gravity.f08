module mod_gravity
  use mod_global_variables, only: dp
  use mod_physics_utils, only: physics_i, physics_zero_func
  use mod_settings, only: settings_t
  use mod_background, only: background_t
  implicit none

  private

  type, public :: gravity_t
    procedure(physics_i), pointer, nopass :: g0
  contains
    procedure :: delete
  end type gravity_t

  public :: new_gravity

contains

  function new_gravity() result(gravity)
    type(gravity_t) :: gravity
    gravity%g0 => physics_zero_func
  end function new_gravity


  pure subroutine delete(this)
    class(gravity_t), intent(inout) :: this
    nullify(this%g0)
  end subroutine delete

end module mod_gravity
