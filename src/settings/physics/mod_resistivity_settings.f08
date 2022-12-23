module mod_resistivity_settings
  use mod_global_variables, only: dp
  implicit none

  private

  type, public :: resistivity_settings_t
    logical, private :: has_resistivity
    real(dp), private :: fixed_eta_value

    logical :: use_fixed_resistivity
    logical :: use_dropoff

  contains

    procedure, public :: enable
    procedure, public :: disable
    procedure, public :: is_enabled
    procedure, public :: set_fixed_resistivity
    procedure, public :: get_fixed_resistivity
    procedure, public :: set_defaults
  end type resistivity_settings_t

  public :: new_resistivity_settings

contains

  pure function new_resistivity_settings() result(resistivity)
    type(resistivity_settings_t) :: resistivity
    call resistivity%set_defaults()
  end function new_resistivity_settings


  pure logical function is_enabled(this)
    class(resistivity_settings_t), intent(in) :: this
    is_enabled = this%has_resistivity
  end function is_enabled


  pure subroutine enable(this)
    class(resistivity_settings_t), intent(inout) :: this
    this%has_resistivity = .true.
  end subroutine enable


  pure subroutine disable(this)
    class(resistivity_settings_t), intent(inout) :: this
    this%has_resistivity = .false.
  end subroutine disable


  pure subroutine set_fixed_resistivity(this, eta)
    class(resistivity_settings_t), intent(inout) :: this
    real(dp), intent(in) :: eta
    this%fixed_eta_value = eta
    this%use_fixed_resistivity = .true.
  end subroutine set_fixed_resistivity


  pure real(dp) function get_fixed_resistivity(this)
    class(resistivity_settings_t), intent(in) :: this
    get_fixed_resistivity = this%fixed_eta_value
  end function get_fixed_resistivity


  pure subroutine set_defaults(this)
    class(resistivity_settings_t), intent(inout) :: this
    this%has_resistivity = .false.
    this%fixed_eta_value = 0.0_dp
    this%use_fixed_resistivity = .false.
    this%use_dropoff = .false.
  end subroutine set_defaults

end module mod_resistivity_settings
