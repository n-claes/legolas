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
  use mod_cooling_curve_names, only: ROSNER
  implicit none

  private

  type(settings_t), pointer :: settings => null()
  type(background_t), pointer :: background => null()

  type, public :: cooling_t
    procedure(real(dp)), pointer, nopass :: lambdaT
    procedure(real(dp)), pointer, nopass :: dlambdadT
    procedure(real(dp)), pointer, nopass :: L0
    procedure(real(dp)), pointer, nopass :: dLdT
    procedure(real(dp)), pointer, nopass :: dLdrho
    logical, private :: is_initialised

  contains
    procedure, public :: initialise
    procedure, public :: delete
  end type cooling_t

  public :: new_cooling

contains

  function new_cooling(settings_tgt, background_tgt) result(cooling)
    type(settings_t), target, intent(in) :: settings_tgt
    type(background_t), target, intent(in) :: background_tgt
    type(cooling_t) :: cooling
    settings => settings_tgt
    background => background_tgt
    cooling%is_initialised = .false.
    cooling%lambdaT => get_lambdaT
    cooling%dlambdadT => get_dlambdadT
    cooling%L0 => get_L0
    cooling%dLdT => get_dLdT
    cooling%dLdrho => get_dLdrho
  end function new_cooling


  subroutine initialise(this)
    use mod_cooling_curves, only: is_valid_cooling_curve, interpolate_cooling_curves

    class(cooling_t), intent(inout) :: this
    character(:), allocatable :: curve
    logical :: use_interpolated_curve

    if (this%is_initialised) return

    curve = settings%physics%cooling%get_cooling_curve()
    if (.not. is_valid_cooling_curve(curve)) return
    use_interpolated_curve = (curve /= ROSNER)
    deallocate(curve)

    if (use_interpolated_curve) then
      call interpolate_cooling_curves(settings)
      this%is_initialised = .true.
    end if
  end subroutine initialise


  real(dp) function get_L0(x)
    real(dp), intent(in) :: x

    get_L0 = 0.0_dp
  end function get_L0


  real(dp) function get_dLdT(x)
    real(dp), intent(in) :: x

    get_dLdT = 0.0_dp
    if (.not. settings%physics%cooling%is_enabled()) return

    get_dLdT = background%density%rho0(x) * get_dlambdadT(x)
  end function get_dLdT


  real(dp) function get_dLdrho(x)
    real(dp), intent(in) :: x

    get_dLdrho = 0.0_dp
    if (.not. settings%physics%cooling%is_enabled()) return

    get_dLdrho = get_lambdaT(x)
  end function get_dLdrho


  real(dp) function get_lambdaT(x)
    use mod_cooling_curves, only: get_rosner_lambdaT, get_interpolated_lambdaT
    real(dp), intent(in) :: x

    if (settings%physics%cooling%get_cooling_curve() == ROSNER) then
      get_lambdaT = get_rosner_lambdaT(x, settings, background)
    else
      get_lambdaT = get_interpolated_lambdaT(x, settings, background)
    end if
  end function get_lambdaT


  real(dp) function get_dlambdadT(x)
    use mod_cooling_curves, only: get_rosner_dlambdadT, get_interpolated_dlambdadT
    real(dp), intent(in) :: x

    if (settings%physics%cooling%get_cooling_curve() == ROSNER) then
      get_dlambdadT = get_rosner_dlambdadT(x, settings, background)
    else
      get_dlambdadT = get_interpolated_dlambdadT(x, settings, background)
    end if
  end function get_dlambdadT


  subroutine delete(this)
    use mod_cooling_curves, only: deallocate_cooling_curves

    class(cooling_t), intent(inout) :: this
    nullify(settings)
    nullify(background)
    nullify(this%lambdaT)
    nullify(this%dlambdadT)
    nullify(this%L0)
    nullify(this%dLdT)
    nullify(this%dLdrho)
    this%is_initialised = .false.
    call deallocate_cooling_curves()
  end subroutine delete

end module mod_radiative_cooling
