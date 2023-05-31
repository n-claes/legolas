module mod_physics_settings
  use mod_global_variables, only: dp
  use mod_check_values, only: is_zero
  use mod_flow_settings, only: flow_settings_t, new_flow_settings
  use mod_cooling_settings, only: cooling_settings_t, new_cooling_settings
  use mod_heating_settings, only: heating_settings_t, new_heating_settings
  use mod_gravity_settings, only: gravity_settings_t, new_gravity_settings
  use mod_resistivity_settings, only: resistivity_settings_t, new_resistivity_settings
  use mod_viscosity_settings, only: viscosity_settings_t, new_viscosity_settings
  use mod_conduction_settings, only: conduction_settings_t, new_conduction_settings
  use mod_hall_settings, only: hall_settings_t, new_hall_settings
  implicit none

  private

  type, public :: physics_settings_t
    real(dp), private :: gamma
    logical :: is_incompressible
    real(dp) :: dropoff_edge_dist
    real(dp) :: dropoff_width
    type(flow_settings_t) :: flow
    type(cooling_settings_t) :: cooling
    type(heating_settings_t) :: heating
    type(gravity_settings_t) :: gravity
    type(resistivity_settings_t) :: resistivity
    type(viscosity_settings_t) :: viscosity
    type(conduction_settings_t) :: conduction
    type(hall_settings_t) :: hall

  contains

    procedure, public :: set_gamma
    procedure, public :: get_gamma
    procedure, public :: get_gamma_1
    procedure, public :: set_incompressible

    procedure, public :: enable_flow
    procedure, public :: enable_cooling
    procedure, public :: enable_heating
    procedure, public :: enable_gravity
    procedure, public :: enable_resistivity
    procedure, public :: enable_viscosity
    procedure, public :: enable_parallel_conduction
    procedure, public :: enable_perpendicular_conduction
    procedure, public :: enable_hall
  end type physics_settings_t

  public :: new_physics_settings

contains

  pure function new_physics_settings() result(physics_settings)
    type(physics_settings_t) :: physics_settings

    call physics_settings%set_gamma(5.0_dp / 3.0_dp)
    physics_settings%is_incompressible = .false.
    physics_settings%dropoff_edge_dist = 0.0_dp
    physics_settings%dropoff_width = 0.0_dp

    physics_settings%flow = new_flow_settings()
    physics_settings%cooling = new_cooling_settings()
    physics_settings%heating = new_heating_settings()
    physics_settings%gravity = new_gravity_settings()
    physics_settings%resistivity = new_resistivity_settings()
    physics_settings%viscosity = new_viscosity_settings()
    physics_settings%conduction = new_conduction_settings()
    physics_settings%hall = new_hall_settings()
  end function new_physics_settings


  pure subroutine set_gamma(this, gamma)
    class(physics_settings_t), intent(inout) :: this
    real(dp), intent(in) :: gamma
    this%gamma = gamma
  end subroutine set_gamma


  pure real(dp) function get_gamma(this)
    class(physics_settings_t), intent(in) :: this
    get_gamma = this%gamma
  end function get_gamma


  pure real(dp) function get_gamma_1(this)
    class(physics_settings_t), intent(in) :: this
    get_gamma_1 = this%get_gamma() - 1.0_dp
  end function get_gamma_1


  pure subroutine set_incompressible(this)
    class(physics_settings_t), intent(inout) :: this
    this%is_incompressible = .true.
    call this%set_gamma(1.0e12_dp)
  end subroutine set_incompressible


  pure subroutine enable_flow(this)
    class(physics_settings_t), intent(inout) :: this
    call this%flow%enable()
  end subroutine enable_flow


  pure subroutine enable_cooling(this, cooling_curve, interpolation_points)
    class(physics_settings_t), intent(inout) :: this
    character(len=*), intent(in), optional :: cooling_curve
    integer, intent(in), optional :: interpolation_points

    if (present(interpolation_points)) then
      call this%cooling%set_interpolation_points(interpolation_points)
    end if
    if (present(cooling_curve)) call this%cooling%set_cooling_curve(cooling_curve)
    call this%cooling%enable()
  end subroutine enable_cooling


  pure subroutine enable_heating(this, force_thermal_balance)
    class(physics_settings_t), intent(inout) :: this
    logical, intent(in), optional :: force_thermal_balance

    if (present(force_thermal_balance)) then
      this%heating%force_thermal_balance = force_thermal_balance
    end if
    call this%heating%enable()
  end subroutine enable_heating


  pure subroutine enable_gravity(this)
    class(physics_settings_t), intent(inout) :: this
    call this%gravity%enable()
  end subroutine enable_gravity


  pure subroutine enable_resistivity(this, fixed_resistivity_value)
    class(physics_settings_t), intent(inout) :: this
    real(dp), intent(in), optional :: fixed_resistivity_value
    real(dp) :: fixed_eta

    call this%resistivity%enable()

    fixed_eta = 0.0_dp
    if (present(fixed_resistivity_value)) fixed_eta = fixed_resistivity_value
    if (.not. is_zero(fixed_eta)) then
      call this%resistivity%set_fixed_resistivity(fixed_eta)
    end if
  end subroutine enable_resistivity


  pure subroutine enable_viscosity(this, viscosity_value, viscous_heating)
    class(physics_settings_t), intent(inout) :: this
    real(dp), intent(in) :: viscosity_value
    logical, intent(in), optional :: viscous_heating
    logical :: heating

    call this%viscosity%set_viscosity_value(viscosity_value)

    heating = this%viscosity%has_viscous_heating()
    if (present(viscous_heating)) heating = viscous_heating
    if (heating) call this%viscosity%enable_viscous_heating()
  end subroutine enable_viscosity


  pure subroutine enable_parallel_conduction(this, fixed_tc_para_value)
    class(physics_settings_t), intent(inout) :: this
    real(dp), intent(in), optional :: fixed_tc_para_value
    real(dp) :: fixed_tc_para

    call this%conduction%enable_para_conduction()

    fixed_tc_para = this%conduction%get_fixed_tc_para()
    if (present(fixed_tc_para_value)) fixed_tc_para = fixed_tc_para_value
    if (.not. is_zero(fixed_tc_para)) then
      call this%conduction%set_fixed_tc_para(fixed_tc_para)
    end if
  end subroutine enable_parallel_conduction


  pure subroutine enable_perpendicular_conduction(this, fixed_tc_perp_value)
    class(physics_settings_t), intent(inout) :: this
    real(dp), intent(in), optional :: fixed_tc_perp_value
    real(dp) :: fixed_tc_perp

    call this%conduction%enable_perp_conduction()

    fixed_tc_perp = this%conduction%get_fixed_tc_perp()
    if (present(fixed_tc_perp_value)) fixed_tc_perp = fixed_tc_perp_value
    if (.not. is_zero(fixed_tc_perp)) then
      call this%conduction%set_fixed_tc_perp(fixed_tc_perp)
    end if
  end subroutine enable_perpendicular_conduction


  pure subroutine enable_hall(this, electron_inertia, electron_fraction)
    class(physics_settings_t), intent(inout) :: this
    logical, intent(in), optional :: electron_inertia
    real(dp), intent(in), optional :: electron_fraction
    logical :: inertia

    call this%hall%enable()

    inertia = this%hall%has_electron_inertia()
    if (present(electron_inertia)) inertia = electron_inertia
    if (inertia) call this%hall%enable_electron_inertia()

    if (present(electron_fraction)) then
      call this%hall%set_electron_fraction(electron_fraction)
    end if
  end subroutine enable_hall

end module mod_physics_settings
