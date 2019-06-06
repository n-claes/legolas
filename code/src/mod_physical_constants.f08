module mod_physical_constants
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter        :: dp = real64

  !! Normalisation-related parameters
  !> Scaling factor for length
  real(dp), protected       :: unit_length = 1.0d0
  !> Scaling factor for time
  real(dp), protected       :: unit_time = 1.0d0
  !> Scaling factor for density
  real(dp), protected       :: unit_density = 1.0d0
  !> Scaling factor for velocity
  real(dp), protected       :: unit_velocity = 0.0d0
  !> Scaling factor for temperature
  real(dp), protected       :: unit_temperature = 1.0d0
  !> Scaling factor for pressure
  real(dp), protected       :: unit_pressure = 1.0d0
  !> Scaling factor for magnetic fields
  real(dp), protected       :: unit_magneticfield = 1.0d0
  !> Scaling factor for number density
  real(dp), protected       :: unit_numberdensity = 1.0d0
  !> Scaling factor for luminosity (radiative cooling module)
  real(dp), protected       :: unit_luminosity = 1.0d0

  !! Physical constants
  !> Value for pi
  real(dp), parameter       :: dpi = 3.141592653589793238462643383279
  !> Helium abundance [He]/[H]
  real(dp), parameter       :: He_abundance = 0.1d0
  !> Proton mass in g (cgs)
  real(dp), parameter       :: mp_cgs = 1.672621777d-24
  !> Proton mass in kg (SI)
  real(dp), parameter       :: mp_si  = 1.672621777d-27
  !> Hydrogen mass in g (cgs)
  real(dp), parameter       :: mH_cgs = 1.6733d-24
  !> Hydrogen mass in kg (SI)
  real(dp), parameter       :: mH_si  = 1.6733d-27
  !> Boltzmann constant in erg/K (cgs)
  real(dp), parameter       :: kB_cgs = 1.3806488d-16
  !> Boltzmann constant in J/K (SI)
  real(dp), parameter       :: kB_si  = 1.3806488d-23
  !> Magnetic constant in H/m (SI)
  real(dp), parameter       :: mu0_si = 1.2566370614d-6

contains

  subroutine set_unit_length(new_unit_length)
    real(dp), intent(in) :: new_unit_length

    unit_length = new_unit_length

  end subroutine set_unit_length


  subroutine set_unit_numberdensity(new_unit_numberdensity)
    real(dp), intent(in) :: new_unit_numberdensity

    unit_numberdensity = new_unit_numberdensity

  end subroutine set_unit_numberdensity


  subroutine set_unit_temperature(new_unit_temperature)
    real(dp), intent(in) :: new_unit_temperature

    unit_temperature = new_unit_temperature

  end subroutine set_unit_temperature


  subroutine set_unit_velocity(new_unit_velocity)
    real(dp), intent(in) :: new_unit_velocity

    unit_velocity = new_unit_velocity

  end subroutine set_unit_velocity


  subroutine set_normalisations(cgs_units)
    logical, intent(in)       :: cgs_units
    real(dp)              :: mp, kB, mu0

    if (cgs_units) then
      mp  = mp_cgs
      kB  = kB_cgs
      mu0 = 4.0d0*dpi
    else
      mp  = mp_si
      kB  = kB_si
      mu0 = mu0_si
    end if

    unit_density  = (1.0d0 + 4.0d0*He_abundance) * mp * unit_numberdensity

    ! if unit_velocity = 0 (L, nb and T defined) -- floating point comparison
    if (abs(unit_velocity - 0.0d0) < 1.0d-9) then
      unit_pressure    = (2.0d0 + 3.0d0*He_abundance) * &
                         unit_numberdensity * kB * unit_temperature
      unit_velocity    = sqrt(unit_pressure / unit_density)
    ! unit velocity != 0 (L, nb and v defined)
    else
      unit_pressure    = unit_density * unit_velocity**2
      unit_temperature = unit_pressure / ((2.0d0 + 3.0d0*He_abundance) * &
                                          unit_numberdensity * kB)
    end if

    unit_magneticfield = sqrt(mu0 * unit_pressure)
    unit_time          = unit_length / unit_velocity

    unit_luminosity    = unit_pressure / ((1.0d0 + 2.0d0*He_abundance) * &
                                          (unit_numberdensity**2 * unit_time))

  end subroutine set_normalisations





end module mod_physical_constants
