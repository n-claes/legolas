module mod_physical_constants
  use, intrinsic :: iso_fortran_env
  implicit none

  ! note: because dp is defined in mod_global_variables, use real64 from
  ! iso_fortran_env directly in this module. It's best to keep this module
  ! without any dependencies and to keep sp, dp and qp in mod_global_variables
  ! (to prevent having to import both modules almost everywhere).

  !! Normalisation-related parameters
  !> Scaling factor for length
  real(real64), protected       :: unit_length = 1.0d0
  !> Scaling factor for time
  real(real64), protected       :: unit_time = 1.0d0
  !> Scaling factor for density
  real(real64), protected       :: unit_density = 1.0d0
  !> Scaling factor for velocity
  real(real64), protected       :: unit_velocity = 0.0d0
  !> Scaling factor for temperature
  real(real64), protected       :: unit_temperature = 1.0d0
  !> Scaling factor for pressure
  real(real64), protected       :: unit_pressure = 1.0d0
  !> Scaling factor for magnetic fields
  real(real64), protected       :: unit_magneticfield = 1.0d0
  !> Scaling factor for number density
  real(real64), protected       :: unit_numberdensity = 1.0d0
  !> Scaling factor for luminosity (radiative cooling module)
  real(real64), protected       :: unit_luminosity = 1.0d0
  !> Scaling factor for themal conduction
  real(real64), protected       :: unit_conduction = 1.0d0

  !! Physical constants
  !> Value for pi
  real(real64), parameter       :: dpi = 3.141592653589793238462643383279
  !> Helium abundance [He]/[H]
  real(real64), parameter       :: He_abundance = 0.1d0
  !> Proton mass in g (cgs)
  real(real64), parameter       :: mp_cgs = 1.672621777d-24
  !> Proton mass in kg (SI)
  real(real64), parameter       :: mp_si  = 1.672621777d-27
  !> Hydrogen mass in g (cgs)
  real(real64), parameter       :: mH_cgs = 1.6733d-24
  !> Hydrogen mass in kg (SI)
  real(real64), parameter       :: mH_si  = 1.6733d-27
  !> Boltzmann constant in erg/K (cgs)
  real(real64), parameter       :: kB_cgs = 1.3806488d-16
  !> Boltzmann constant in J/K (SI)
  real(real64), parameter       :: kB_si  = 1.3806488d-23
  !> Magnetic constant in H/m (SI)
  real(real64), parameter       :: mu0_si = 1.2566370614d-6

contains

  subroutine set_unit_length(new_unit_length)
    real(real64), intent(in) :: new_unit_length

    unit_length = new_unit_length

  end subroutine set_unit_length


  subroutine set_unit_numberdensity(new_unit_numberdensity)
    real(real64), intent(in) :: new_unit_numberdensity

    unit_numberdensity = new_unit_numberdensity

  end subroutine set_unit_numberdensity


  subroutine set_unit_temperature(new_unit_temperature)
    real(real64), intent(in) :: new_unit_temperature

    unit_temperature = new_unit_temperature

  end subroutine set_unit_temperature


  subroutine set_unit_velocity(new_unit_velocity)
    real(real64), intent(in) :: new_unit_velocity

    unit_velocity = new_unit_velocity

  end subroutine set_unit_velocity


  subroutine set_normalisations(cgs_units)
    logical, intent(in)       :: cgs_units
    real(real64)              :: mp, kB, mu0

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

    unit_conduction    = unit_density * unit_length * unit_velocity**3 / &
                         unit_temperature

  end subroutine set_normalisations





end module mod_physical_constants
