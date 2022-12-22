module mod_viscosity_settings
  use mod_global_variables, only: dp
  implicit none

  private

  type, public :: viscosity_settings_t
    logical, private :: has_viscosity
    logical, private :: viscous_heating
    real(dp), private :: fixed_viscosity_value

  contains

    procedure, public :: enable
    procedure, public :: enable_viscous_heating
    procedure, public :: disable
    procedure, public :: is_enabled
    procedure, public :: has_viscous_heating
    procedure, public :: set_fixed_viscosity
    procedure, public :: get_fixed_viscosity
  end type viscosity_settings_t

  public :: new_viscosity_settings

contains

  pure function new_viscosity_settings() result(viscosity_settings)
    type(viscosity_settings_t) :: viscosity_settings
    viscosity_settings%has_viscosity = .false.
    viscosity_settings%viscous_heating = .false.
    viscosity_settings%fixed_viscosity_value = 0.0_dp
  end function new_viscosity_settings


  pure logical function is_enabled(this)
    class(viscosity_settings_t), intent(in) :: this
    is_enabled = this%has_viscosity
  end function is_enabled


  pure subroutine enable(this)
    class(viscosity_settings_t), intent(inout) :: this
    this%has_viscosity = .true.
  end subroutine enable


  pure subroutine enable_viscous_heating(this)
    class(viscosity_settings_t), intent(inout) :: this
    this%viscous_heating = .true.
  end subroutine enable_viscous_heating


  pure subroutine disable(this)
    class(viscosity_settings_t), intent(inout) :: this
    this%has_viscosity = .false.
    this%viscous_heating = .false.
  end subroutine disable


  pure logical function has_viscous_heating(this)
    class(viscosity_settings_t), intent(in) :: this
    has_viscous_heating = this%viscous_heating
  end function has_viscous_heating


  pure subroutine set_fixed_viscosity(this, viscosity_value)
    class(viscosity_settings_t), intent(inout) :: this
    real(dp), intent(in) :: viscosity_value
    this%fixed_viscosity_value = viscosity_value
  end subroutine set_fixed_viscosity


  pure real(dp) function get_fixed_viscosity(this)
    class(viscosity_settings_t), intent(in) :: this
    get_fixed_viscosity = this%fixed_viscosity_value
  end function get_fixed_viscosity

end module mod_viscosity_settings
