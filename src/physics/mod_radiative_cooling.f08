! =============================================================================
!> Module containing radiative cooling-related routines.
!! This module is responsible for initialising the radiative cooling
!! variables and a correct handling of the cooling curves.
!! If an interpolated cooling curve is selected this module calls the
!! interpolation module to create one.
module mod_radiative_cooling
  use mod_global_variables, only: dp
  use mod_logging, only: logger
  use mod_settings, only: settings_t
  use mod_background, only: background_t
  use mod_physics_utils, only: physics_i
  use mod_cooling_curves, only: interpolate_cooling_curves, get_lambdaT, get_dlambdadT
  implicit none

  private

  type, public :: cooling_t
    procedure(physics_i), pointer, nopass :: lambdaT
    procedure(physics_i), pointer, nopass :: dlambdadT
    procedure(physics_i), pointer, nopass :: L0
    procedure(physics_i), pointer, nopass :: dLdT
    procedure(physics_i), pointer, nopass :: dLdrho
    logical, private :: is_initialised

  contains
    procedure, public :: initialise
    procedure, public :: delete
  end type cooling_t

  public :: new_cooling

contains

  function new_cooling() result(cooling)
    type(cooling_t) :: cooling
    cooling%is_initialised = .false.
    cooling%lambdaT => get_lambdaT
    cooling%dlambdadT => get_dlambdadT
    cooling%L0 => get_L0
    cooling%dLdT => get_dLdT
    cooling%dLdrho => get_dLdrho
  end function new_cooling


  real(dp) function zero(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    zero = 0.0_dp
  end function zero


  subroutine initialise(this, settings)
    class(cooling_t), intent(inout) :: this
    type(settings_t), intent(in) :: settings
    logical :: use_interpolated_curve

    if (this%is_initialised) return

    use_interpolated_curve = (settings%physics%cooling%get_cooling_curve() /= "rosner")
    if (use_interpolated_curve) then
      call interpolate_cooling_curves(settings)
      this%is_initialised = .true.
    end if
  end subroutine initialise


  real(dp) function get_L0(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background

    get_L0 = 0.0_dp
  end function get_L0


  real(dp) function get_dLdT(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background

    get_dLdT = 0.0_dp
    if (.not. settings%physics%cooling%is_enabled()) return

    get_dLdT = background%density%rho0(x) * get_dlambdadT(x, settings, background)
  end function get_dLdT


  real(dp) function get_dLdrho(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background

    get_dLdrho = 0.0_dp
    if (.not. settings%physics%cooling%is_enabled()) return

    get_dLdrho = get_lambdaT(x, settings, background)
  end function get_dLdrho


  subroutine delete(this)
    use mod_cooling_curves, only: deallocate_cooling_curves

    class(cooling_t), intent(inout) :: this
    nullify(this%lambdaT)
    nullify(this%dlambdadT)
    nullify(this%L0)
    nullify(this%dLdT)
    nullify(this%dLdrho)
    this%is_initialised = .false.
    call deallocate_cooling_curves()
  end subroutine delete

end module mod_radiative_cooling
