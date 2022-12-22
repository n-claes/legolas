module mod_physics_settings
  use mod_global_variables, only: dp
  use mod_flow_settings, only: flow_settings_t, new_flow_settings
  use mod_cooling_settings, only: cooling_settings_t, new_cooling_settings
  use mod_gravity_settings, only: gravity_settings_t, new_gravity_settings
  implicit none

  private

  type, public :: physics_t
    real(dp), private :: gamma
    logical :: is_incompressible
    type(flow_settings_t) :: flow
    type(cooling_settings_t) :: cooling
    type(gravity_settings_t) :: gravity

  contains

    procedure, public :: set_gamma
    procedure, public :: get_gamma
    procedure, public :: get_gamma_1
    procedure, public :: set_incompressible
    procedure, public :: set_defaults

    procedure, public :: enable_flow
    procedure, public :: enable_cooling
    procedure, public :: enable_gravity
  end type physics_t

  public :: new_physics_settings

contains

  pure function new_physics_settings() result(physics)
    type(physics_t) :: physics

    physics%is_incompressible = .false.
  end function new_physics_settings


  pure subroutine set_gamma(this, gamma)
    class(physics_t), intent(inout) :: this
    real(dp), intent(in) :: gamma
    this%gamma = gamma
  end subroutine set_gamma


  pure real(dp) function get_gamma(this)
    class(physics_t), intent(in) :: this
    get_gamma = this%gamma
  end function get_gamma


  pure real(dp) function get_gamma_1(this)
    class(physics_t), intent(in) :: this
    get_gamma_1 = this%get_gamma() - 1.0_dp
  end function get_gamma_1


  pure subroutine set_incompressible(this)
    class(physics_t), intent(inout) :: this
    this%is_incompressible = .true.
    call this%set_gamma(1.0e12_dp)
  end subroutine set_incompressible


  pure subroutine set_defaults(this)
    class(physics_t), intent(inout) :: this
    call this%set_gamma(5.0_dp / 3.0_dp)
    this%is_incompressible = .false.
  end subroutine set_defaults


  pure subroutine enable_flow(this)
    class(physics_t), intent(inout) :: this
    call this%flow%enable()
  end subroutine enable_flow


  pure subroutine enable_cooling(this, cooling_curve, interpolation_points)
    class(physics_t), intent(inout) :: this
    character(len=*), intent(in) :: cooling_curve
    integer, intent(in), optional :: interpolation_points

    if (present(interpolation_points)) then
      call this%cooling%set_interpolation_points(interpolation_points)
    end if
    call this%cooling%set_cooling_curve(cooling_curve)
    call this%cooling%enable()
  end subroutine enable_cooling


  pure subroutine enable_gravity(this)
    class(physics_t), intent(inout) :: this
    call this%gravity%enable()
  end subroutine enable_gravity

end module mod_physics_settings
