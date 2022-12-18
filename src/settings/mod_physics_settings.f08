module mod_physics_settings
  use mod_global_variables, only: dp
  implicit none

  private

  type, public :: physics_t
    real(dp), private :: gamma
    logical :: is_incompressible

  contains

    procedure, public :: set_gamma
    procedure, public :: get_gamma
    procedure, public :: get_gamma_1
    procedure, public :: set_incompressible
    procedure, public :: set_defaults
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

end module mod_physics_settings