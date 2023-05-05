module mod_resistivity_settings
  use mod_global_variables, only: dp
  implicit none

  private

  type, public :: resistivity_settings_t
    logical, private :: has_resistivity
    real(dp), private :: fixed_resistivity_value
    logical, private :: fixed_resistivity
    logical :: use_dropoff

  contains

    procedure, public :: enable
    procedure, public :: disable
    procedure, public :: is_enabled
    procedure, public :: set_fixed_resistivity
    procedure, public :: get_fixed_resistivity
    procedure, public :: has_fixed_resistivity
  end type resistivity_settings_t

  public :: new_resistivity_settings

contains

  pure function new_resistivity_settings() result(resistivity)
    type(resistivity_settings_t) :: resistivity

    resistivity%has_resistivity = .false.
    resistivity%fixed_resistivity_value = 0.0_dp
    resistivity%fixed_resistivity = .false.
    resistivity%use_dropoff = .false.
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
    this%fixed_resistivity_value = eta
    this%fixed_resistivity = .true.
    this%has_resistivity = .true.
  end subroutine set_fixed_resistivity


  pure real(dp) function get_fixed_resistivity(this)
    class(resistivity_settings_t), intent(in) :: this
    get_fixed_resistivity = this%fixed_resistivity_value
  end function get_fixed_resistivity


  pure logical function has_fixed_resistivity(this)
    class(resistivity_settings_t), intent(in) :: this
    has_fixed_resistivity = this%fixed_resistivity
  end function has_fixed_resistivity

end module mod_resistivity_settings
