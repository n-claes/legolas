module mod_bg_velocity
  use mod_global_variables, only: dp
  implicit none

  private

  type, public :: bg_velocity_t
    procedure(real(dp)), pointer, nopass :: v01
    procedure(real(dp)), pointer, nopass :: dv01
    procedure(real(dp)), pointer, nopass :: ddv01
    procedure(real(dp)), pointer, nopass :: v02
    procedure(real(dp)), pointer, nopass :: dv02
    procedure(real(dp)), pointer, nopass :: ddv02
    procedure(real(dp)), pointer, nopass :: v03
    procedure(real(dp)), pointer, nopass :: dv03
    procedure(real(dp)), pointer, nopass :: ddv03
  contains
    procedure, public :: get_v0
    procedure, public :: delete
  end type bg_velocity_t

  public :: new_bg_velocity

contains

  function new_bg_velocity(default_func) result(bg_velocity)
    procedure(real(dp)) :: default_func
    type(bg_velocity_t) :: bg_velocity
    bg_velocity%v01 => default_func
    bg_velocity%dv01 => default_func
    bg_velocity%ddv01 => default_func
    bg_velocity%v02 => default_func
    bg_velocity%dv02 => default_func
    bg_velocity%ddv02 => default_func
    bg_velocity%v03 => default_func
    bg_velocity%dv03 => default_func
    bg_velocity%ddv03 => default_func
  end function new_bg_velocity


  real(dp) function get_v0(this, x)
    class(bg_velocity_t), intent(in) :: this
    real(dp), intent(in) :: x
    get_v0 = sqrt(this%v01(x)**2 + this%v02(x)**2 + this%v03(x)**2)
  end function get_v0


  pure subroutine delete(this)
    class(bg_velocity_t), intent(inout) :: this
    nullify(this%v01)
    nullify(this%dv01)
    nullify(this%ddv01)
    nullify(this%v02)
    nullify(this%dv02)
    nullify(this%ddv02)
    nullify(this%v03)
    nullify(this%dv03)
    nullify(this%ddv03)
  end subroutine delete

end module mod_bg_velocity
