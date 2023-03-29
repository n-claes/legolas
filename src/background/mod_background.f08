module mod_background
  use mod_global_variables, only: dp
  use mod_function_utils, only: zero_func
  use mod_bg_density, only: bg_density_t, new_bg_density
  use mod_bg_velocity, only: bg_velocity_t, new_bg_velocity
  use mod_bg_temperature, only: bg_temperature_t, new_bg_temperature
  use mod_bg_magnetic, only: bg_magnetic_t, new_bg_magnetic
  implicit none

  private

  type, public :: background_t
    type(bg_density_t) :: density
    type(bg_velocity_t) :: velocity
    type(bg_temperature_t) :: temperature
    type(bg_magnetic_t) :: magnetic

  contains

    procedure :: set_density_funcs
    procedure :: set_velocity_1_funcs
    procedure :: set_velocity_2_funcs
    procedure :: set_velocity_3_funcs
    procedure :: set_temperature_funcs
    procedure :: set_magnetic_1_funcs
    procedure :: set_magnetic_2_funcs
    procedure :: set_magnetic_3_funcs
    procedure :: delete
  end type background_t

  public :: new_background

contains

  function new_background() result(background)
    type(background_t) :: background

    background%density = new_bg_density(default_func=zero_func)
    background%velocity = new_bg_velocity(default_func=zero_func)
    background%temperature = new_bg_temperature(default_func=zero_func)
    background%magnetic = new_bg_magnetic(default_func=zero_func)
  end function new_background

  subroutine set_density_funcs(this, rho0_func, drho0_func)
    class(background_t), intent(inout) :: this
    procedure(real(dp)) :: rho0_func
    procedure(real(dp)), optional :: drho0_func

    this%density%rho0 => rho0_func
    if (present(drho0_func)) this%density%drho0 => drho0_func
  end subroutine set_density_funcs


  subroutine set_velocity_1_funcs(this, v01_func, dv01_func, ddv01_func)
    class(background_t), intent(inout) :: this
    procedure(real(dp)) :: v01_func
    procedure(real(dp)), optional :: dv01_func
    procedure(real(dp)), optional :: ddv01_func

    this%velocity%v01 => v01_func
    if (present(dv01_func)) this%velocity%dv01 => dv01_func
    if (present(ddv01_func)) this%velocity%ddv01 => ddv01_func
  end subroutine set_velocity_1_funcs


  subroutine set_velocity_2_funcs(this, v02_func, dv02_func, ddv02_func)
    class(background_t), intent(inout) :: this
    procedure(real(dp)) :: v02_func
    procedure(real(dp)), optional :: dv02_func
    procedure(real(dp)), optional :: ddv02_func

    this%velocity%v02 => v02_func
    if (present(dv02_func)) this%velocity%dv02 => dv02_func
    if (present(ddv02_func)) this%velocity%ddv02 => ddv02_func
  end subroutine set_velocity_2_funcs


  subroutine set_velocity_3_funcs(this, v03_func, dv03_func, ddv03_func)
    class(background_t), intent(inout) :: this
    procedure(real(dp)) :: v03_func
    procedure(real(dp)), optional :: dv03_func
    procedure(real(dp)), optional :: ddv03_func

    this%velocity%v03 => v03_func
    if (present(dv03_func)) this%velocity%dv03 => dv03_func
    if (present(ddv03_func)) this%velocity%ddv03 => ddv03_func
  end subroutine set_velocity_3_funcs


  subroutine set_temperature_funcs(this, T0_func, dT0_func, ddT0_func)
    class(background_t), intent(inout) :: this
    procedure(real(dp)) :: T0_func
    procedure(real(dp)), optional :: dT0_func
    procedure(real(dp)), optional :: ddT0_func

    this%temperature%T0 => T0_func
    if (present(dT0_func)) this%temperature%dT0 => dT0_func
    if (present(ddT0_func)) this%temperature%ddT0 => ddT0_func
  end subroutine set_temperature_funcs


  subroutine set_magnetic_1_funcs(this, B01_func)
    class(background_t), intent(inout) :: this
    procedure(real(dp)) :: B01_func

    this%magnetic%B01 => B01_func
  end subroutine set_magnetic_1_funcs


  subroutine set_magnetic_2_funcs(this, B02_func, dB02_func, ddB02_func)
    class(background_t), intent(inout) :: this
    procedure(real(dp)) :: B02_func
    procedure(real(dp)), optional :: dB02_func
    procedure(real(dp)), optional :: ddB02_func

    this%magnetic%B02 => B02_func
    if (present(dB02_func)) this%magnetic%dB02 => dB02_func
    if (present(ddB02_func)) this%magnetic%ddB02 => ddB02_func
  end subroutine set_magnetic_2_funcs


  subroutine set_magnetic_3_funcs(this, B03_func, dB03_func, ddB03_func)
    class(background_t), intent(inout) :: this
    procedure(real(dp)) :: B03_func
    procedure(real(dp)), optional :: dB03_func
    procedure(real(dp)), optional :: ddB03_func

    this%magnetic%B03 => B03_func
    if (present(dB03_func)) this%magnetic%dB03 => dB03_func
    if (present(ddB03_func)) this%magnetic%ddB03 => ddB03_func
  end subroutine set_magnetic_3_funcs


  pure subroutine delete(this)
    class(background_t), intent(inout) :: this
    call this%density%delete()
    call this%temperature%delete()
    call this%magnetic%delete()
    call this%velocity%delete()
  end subroutine delete

end module mod_background
