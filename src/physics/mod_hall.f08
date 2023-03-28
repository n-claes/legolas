! =============================================================================
!> Module containing Hall-related routines.
!! Sets the Hall and electron inertia factors based on normalisations and specified profiles.
module mod_hall
  use mod_global_variables, only: dp
  use mod_settings, only: settings_t
  use mod_background, only: background_t
  use mod_physics_utils, only: physics_i, get_dropoff, get_dropoff_dr
  use mod_physical_constants, only: dpi, mp_cgs, ec_cgs, me_cgs
  implicit none

  private

  type, public :: hall_t
    procedure(physics_i), pointer, nopass :: hallfactor
    procedure(physics_i), pointer, nopass :: inertiafactor

  contains
    procedure, public :: validate_scale_ratio
    procedure, public :: delete
  end type hall_t

  public :: new_hall

contains


  function new_hall() result(hall)
    type(hall_t) :: hall
    hall%hallfactor => get_hallfactor
    hall%inertiafactor => get_inertiafactor
  end function new_hall


  !> Retrieves the normalised Hall factor as described by Porth et al. (2014).
  real(dp) function get_hallfactor(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: unit_velocity, unit_length, unit_magneticfield

    get_hallfactor = 0.0_dp
    if (.not. settings%physics%hall%is_enabled()) return

    unit_velocity = settings%units%get_unit_velocity()
    unit_length = settings%units%get_unit_length()
    unit_magneticfield = settings%units%get_unit_magneticfield()
    get_hallfactor = ( &
      (mp_cgs * unit_velocity) / (ec_cgs * unit_length * unit_magneticfield) &
    )

    if (.not. settings%physics%hall%use_dropoff) return
    get_hallfactor = get_dropoff(x, get_hallfactor, settings, background)
  end function get_hallfactor


  real(dp) function get_inertiafactor(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: unit_velocity, unit_length, unit_magneticfield

    get_inertiafactor = 0.0_dp
    if (.not. settings%physics%hall%has_electron_inertia()) return

    unit_velocity = settings%units%get_unit_velocity()
    unit_length = settings%units%get_unit_length()
    unit_magneticfield = settings%units%get_unit_magneticfield()
    get_inertiafactor = ( &
      mp_cgs * me_cgs * unit_velocity**2 &
      / (ec_cgs * unit_length * unit_magneticfield)**2 &
    )

    if (.not. settings%physics%hall%use_inertia_dropoff) return
    get_inertiafactor = get_dropoff_dr(x, get_inertiafactor, settings, background)
  end function get_inertiafactor


  ! LCOV_EXCL_START
  subroutine validate_scale_ratio(this, grid, settings, background)
    use mod_physics_utils, only: from_physics_function
    use mod_logging, only: logger, str

    class(hall_t), intent(in) :: this
    real(dp), intent(in) :: grid(:)
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: ratio

    if (.not. settings%physics%hall%is_enabled()) return

    ratio = maxval( &
      from_physics_function(this%hallfactor, grid, settings, background) &
    ) / (settings%grid%get_grid_end() - settings%grid%get_grid_start())

    if (ratio > 0.1d0) then
      call logger%warning("large ratio Hall scale / system scale: " // str(ratio))
    end if
  end subroutine validate_scale_ratio
  ! LCOV_EXCL_STOP


  pure subroutine delete(this)
    class(hall_t), intent(inout) :: this
    nullify(this%hallfactor)
    nullify(this%inertiafactor)
  end subroutine delete

end module mod_hall
