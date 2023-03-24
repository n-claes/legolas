module mod_gravity
  use mod_global_variables, only: dp
  use mod_physics_utils, only: physics_i
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
    gravity%g0 => default_g0
  end function new_gravity


  real(dp) function default_g0(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    default_g0 = 0.0_dp
  end function default_g0


  pure subroutine delete(this)
    class(gravity_t), intent(inout) :: this
    nullify(this%g0)
  end subroutine delete

end module mod_gravity
