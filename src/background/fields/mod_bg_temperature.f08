module mod_bg_temperature
  use mod_global_variables, only: dp
  implicit none

  private

  type, public :: bg_temperature_t
    procedure(real(dp)), pointer, nopass :: T0
    procedure(real(dp)), pointer, nopass :: dT0
    procedure(real(dp)), pointer, nopass :: ddT0
  contains
    procedure :: delete
  end type bg_temperature_t

  public :: new_bg_temperature

contains

  function new_bg_temperature(default_func) result(bg_temperature)
    procedure(real(dp)) :: default_func
    type(bg_temperature_t) :: bg_temperature
    bg_temperature%T0 => default_func
    bg_temperature%dT0 => default_func
    bg_temperature%ddT0 => default_func
  end function new_bg_temperature


  pure subroutine delete(this)
    class(bg_temperature_t), intent(inout) :: this
    nullify(this%T0)
    nullify(this%dT0)
    nullify(this%ddT0)
  end subroutine delete

end module mod_bg_temperature
