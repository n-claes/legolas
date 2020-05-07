!
! MODULE: mod_units
!
! DESCRIPTION:
!> Module containing all unit normalisations.
!
module mod_units
  use mod_global_variables, only: dp, cgs_units
  use mod_logging, only: log_message
  implicit none

  private

  !! Normalisation-related parameters
  !> Scaling factor for length
  real(dp), protected   :: unit_length = 1.0d0
  !> Scaling factor for time
  real(dp), protected   :: unit_time = 1.0d0
  !> Scaling factor for density
  real(dp), protected   :: unit_density = 1.0d0
  !> Scaling factor for velocity
  real(dp), protected   :: unit_velocity = 0.0d0
  !> Scaling factor for temperature
  real(dp), protected   :: unit_temperature = 1.0d0
  !> Scaling factor for pressure
  real(dp), protected   :: unit_pressure = 1.0d0
  !> Scaling factor for magnetic fields
  real(dp), protected   :: unit_magneticfield = 1.0d0
  !> Scaling factor for number density
  real(dp), protected   :: unit_numberdensity = 1.0d0
  !> Scaling factor for luminosity (radiative cooling module)
  real(dp), protected   :: unit_luminosity = 1.0d0
  !> Scaling factor for derivative cooling function to temperature
  real(dp), protected   :: unit_dldt = 1.0d0
  !> Scaling factor for themal conduction
  real(dp), protected   :: unit_conduction = 1.0d0
  !> Scaling factor for d(kappa_perp)/d(rho)
  real(dp), protected   :: unit_dtc_drho = 1.0d0
  !> Scaling factor for d(kappa_perp)/d(T)
  real(dp), protected   :: unit_dtc_dT = 1.0d0
  !> Scaling factor for d(kappa_perp)/d(B^2)
  real(dp), protected   :: unit_dtc_dB2 = 1.0d0

  !> Scaling factor for the resistivity eta. This is defined such that the
  !! normalised resistivity at 1 MK equals approximately 0.1.
  real(dp), protected   :: unit_resistivity = 1.0d0
  !> Scaling factor for d(eta)/d(T)
  real(dp), protected   :: unit_deta_dT = 1.0d0

  !> Boolean value to check if normalisations are set
  logical, protected    :: normalisations_are_set = .false.

  public  :: check_if_normalisations_set
  public  :: set_normalisations
  public  :: set_unit_resistivity

  public  :: unit_length, unit_time, unit_density, unit_velocity
  public  :: unit_temperature, unit_pressure, unit_magneticfield
  public  :: unit_numberdensity, unit_luminosity, unit_dldt
  public  :: unit_conduction, unit_dtc_drho, unit_dtc_dT, unit_dtc_dB2
  public  :: unit_resistivity, unit_deta_dT

contains


  !> Checks if normalisations were set. If not, define them using default values, taken
  !! to be unit_temperature = 1 MK, unit_magneticfield = 10 G and unit_length = 1e9 cm.
  subroutine check_if_normalisations_set()
    if (normalisations_are_set) then
      call log_message("normalisations are already set", level='debug')
      return
    else
      cgs_units = .true.
      call set_normalisations(new_unit_temperature=1.0d6, new_unit_magneticfield=10.0d0, new_unit_length=1.0d9)
    end if
  end subroutine check_if_normalisations_set


  !> Returns the Boltzmann constant, proton mass and magnetic permeability either in SI or cgs,
  !! based on the value of use_cgs.
  !! @param[out]  kB  Boltzmann constant
  !! @param[out]  mp  proton mass
  !! @param[out]  mu0 magnetic permeability
  !! @param[out]  Rgas gas constant
  subroutine get_constants(kB, mp, mu0, Rgas)
    use mod_physical_constants, only: kB_si, kB_cgs, mp_si, mp_cgs, mu0_si, mu0_cgs, R_si, R_cgs

    real(dp), intent(out) :: kB, mp, mu0, Rgas

    if (cgs_units) then
      call log_message("getting constants in cgs units", level='debug')
      kB = kB_cgs
      mp = mp_cgs
      mu0 = mu0_cgs
      Rgas = R_cgs
    else
      call log_message("getting constants in SI units", level='debug')
      kB = kB_si
      mp = mp_si
      mu0 = mu0_si
      Rgas = R_si
    end if
  end subroutine get_constants


  !> Defines unit normalisations based on a magnetic field unit, length unit, and a density OR temperature unit.
  !! When calling this routine you should specify the keyword arguments.
  !! Normalisations use a reference value of 2 for plasma beta as is common practice, such that we can write
  !! beta / 2 = mu * p / B**2 = 1. This means that the unit pressure can be written as p = B**2 / 2,
  !! such that temperature or density can be derived using the ideal gas law.
  !! @param[in] new_unit_density        unit density value
  !! @param[in] new_unit_temperature    unit temperature value
  !! @param[in] new_unit_magneticfield  unit magneticfield value
  !! @param[in] new_unit_length         unit length value
  subroutine set_normalisations(new_unit_density, new_unit_temperature, new_unit_magneticfield, new_unit_length)
    real(dp), intent(in), optional  :: new_unit_density, new_unit_temperature
    real(dp), intent(in)            :: new_unit_magneticfield, new_unit_length
    real(dp)  :: kB, mp, mu0, Rgas

    call get_constants(kB, mp, mu0, Rgas)

    if (present(new_unit_density) .and. present(new_unit_temperature)) then
      call log_message("unit density and unit temperature can not both be set.", level='error')
    end if

    unit_magneticfield = new_unit_magneticfield
    unit_length = new_unit_length
    unit_pressure = unit_magneticfield**2 / mu0

    if (present(new_unit_density)) then
      unit_density = new_unit_density
      unit_temperature = unit_pressure / (Rgas * unit_density)
    else if (present(new_unit_temperature)) then
      unit_temperature = new_unit_temperature
      unit_density = unit_pressure / (Rgas * unit_temperature)
    else
      call log_message("no unit density or unit temperature specified.", level='error')
    end if

    unit_numberdensity = unit_density / mp
    unit_velocity = unit_magneticfield / sqrt(mu0 * unit_density)
    unit_time = unit_length / unit_velocity
    unit_luminosity = unit_pressure / (unit_time * unit_numberdensity**2)
    unit_dldt = unit_luminosity / unit_temperature

    unit_conduction    = unit_density * unit_length * unit_velocity**3 / unit_temperature
    unit_dtc_drho      = unit_conduction / unit_density
    unit_dtc_dT        = unit_conduction / unit_temperature
    unit_dtc_dB2       = unit_conduction / (unit_magneticfield**2)

    normalisations_are_set = .true.

  end subroutine


  !> Sets the unit resistivity based on the given value. This is done automatically in the
  !! resistivity module, and is set as such that the normalised resistivity at 1 MK equals 0.1.
  !! @param[in] new_unit_resistivity  resistivity
  subroutine set_unit_resistivity(new_unit_resistivity)
    real(dp), intent(in)  :: new_unit_resistivity

    unit_resistivity = new_unit_resistivity
    unit_deta_dT = unit_resistivity / unit_temperature
  end subroutine set_unit_resistivity

end module mod_units
