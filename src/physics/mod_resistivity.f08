! =============================================================================
!> Module containing resistivity-related routines, calculates
!! and sets the resistivity values based on the equilibrium configuration.
module mod_resistivity
  use mod_global_variables, only: dp
  use mod_physical_constants, only: dpi, Z_ion, coulomb_log, ec_cgs, me_cgs, kB_cgs
  use mod_settings, only: settings_t
  use mod_background, only: background_t
  use mod_physics_utils, only: get_dropoff, get_dropoff_dr
  implicit none

  private

  type(settings_t), pointer :: settings => null()
  type(background_t), pointer :: background => null()

  type, public :: resistivity_t
    procedure(real(dp)), pointer, nopass :: eta
    procedure(real(dp)), pointer, nopass :: detadT
    procedure(real(dp)), pointer, nopass :: detadr
  contains
    procedure, public :: delete
  end type resistivity_t

  public :: new_resistivity

contains


  function new_resistivity(settings_tgt, background_tgt) result(resistivity)
    type(settings_t), target, intent(in) :: settings_tgt
    type(background_t), target, intent(in) :: background_tgt
    type(resistivity_t) :: resistivity

    settings => settings_tgt
    background => background_tgt

    resistivity%eta => spitzer_eta
    resistivity%detadT => spitzer_detadT
    resistivity%detadr => spitzer_detadr
  end function new_resistivity


  !> Default profile for the resistivity $\eta$. Returns one of the following:
  !! - The Spitzer resistivity based on the equilibrium temperature profile
  !! - a fixed resistivity value
  !! - zero (if resistivity is disabled)
  !! Values are normalised on return.
  real(dp) function spitzer_eta(x)
    real(dp), intent(in) :: x
    real(dp) :: unit_temperature, unit_resistivity
    real(dp) :: T0

    spitzer_eta = 0.0_dp
    if (.not. settings%physics%resistivity%is_enabled()) return

    if (settings%physics%resistivity%has_fixed_resistivity()) then
      spitzer_eta = settings%physics%resistivity%get_fixed_resistivity()
      if (settings%physics%resistivity%use_dropoff) then
        spitzer_eta = get_dropoff(x, spitzer_eta, settings)
      end if
      return
    end if

    unit_temperature = settings%units%get_unit_temperature()
    unit_resistivity = settings%units%get_unit_resistivity()
    T0 = background%temperature%T0(x) * unit_temperature
    spitzer_eta = ( &
      (4.0_dp / 3.0_dp) &
      * sqrt(2.0_dp * dpi) &
      * Z_ion &
      * ec_cgs**2 &
      * sqrt(me_cgs) &
      * coulomb_log &
      / (kB_cgs * T0)**(3.0_dp / 2.0_dp) &
    ) / unit_resistivity
  end function spitzer_eta


  !> Default profile for the derivative of the resistivity $\eta$ with respect to
  !! the temperature $T$. Returns one of the following:
  !! - zero (if resistivity is disabled or for fixed resistivity)
  !! - derivative of Spitzer resistivity with respect to T
  !! Values are normalised on return.
  real(dp) function spitzer_detadT(x)
    real(dp), intent(in) :: x
    real(dp) :: unit_temperature, unit_detadT
    real(dp) :: T0

    spitzer_detadT = 0.0_dp
    if (.not. settings%physics%resistivity%is_enabled()) return
    if (settings%physics%resistivity%has_fixed_resistivity()) return

    unit_temperature = settings%units%get_unit_temperature()
    unit_detadT = settings%units%get_unit_resistivity() / unit_temperature
    T0 = background%temperature%T0(x) * unit_temperature
    spitzer_detadT = ( &
      -2.0_dp * sqrt(2.0_dp * dpi) * Z_ion * ec_cgs**2 * sqrt(me_cgs) * coulomb_log &
      / (kB_cgs**(3.0_dp / 2.0_dp) * T0**(5.0_dp / 2.0_dp)) &
    ) / unit_detadT
  end function spitzer_detadT


  !> Default profile for the derivative of the resistivity $\eta$ with respect to
  !! position. Returns zero unless a dropoff profile was chosen.
  real(dp) function spitzer_detadr(x)
    real(dp), intent(in) :: x

    spitzer_detadr = 0.0_dp
    if (.not. settings%physics%resistivity%is_enabled()) return

    if (settings%physics%resistivity%has_fixed_resistivity()) then
      if (settings%physics%resistivity%use_dropoff) then
        spitzer_detadr = get_dropoff_dr( &
          x, settings%physics%resistivity%get_fixed_resistivity(), settings &
        )
      end if
      return
    end if
  end function spitzer_detadr


  subroutine delete(this)
    class(resistivity_t), intent(inout) :: this
    nullify(settings)
    nullify(background)
    nullify(this%eta)
    nullify(this%detadT)
    nullify(this%detadr)
  end subroutine delete

end module mod_resistivity
