!
! MODULE: mod_physical_constants
!
!> @author
!> Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Module containing all physical constants and unit normalisations.
!
module mod_physical_constants
  use, intrinsic :: iso_fortran_env
  implicit none

  public

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
  !> Scaling factor for derivative cooling function to temperature
  real(real64), protected       :: unit_dldt = 1.0d0
  !> Scaling factor for themal conduction
  real(real64), protected       :: unit_conduction = 1.0d0
  !> Scaling factor for d(kappa_perp)/d(rho)
  real(real64), protected       :: unit_dtc_drho = 1.0d0
  !> Scaling factor for d(kappa_perp)/d(T)
  real(real64), protected       :: unit_dtc_dT = 1.0d0
  !> Scaling factor for d(kappa_perp)/d(B^2)
  real(real64), protected       :: unit_dtc_dB2 = 1.0d0
  !> Scaling factor for the resistivity eta. This is defined such that the
  !! normalised resistivity at 1 MK equals approximately 0.1.
  real(real64), protected       :: unit_resistivity = 1.0d0

  !! Physical constants
  !> Value for pi
  real(real64), parameter       :: dpi = 3.141592653589793238462643383279
  !> Helium abundance [He]/[H]
  real(real64), parameter       :: He_abundance = 0.1d0
  !> Coulomb logarithm
  real(real64), parameter       :: coulomb_log = 22
  !> Proton mass in g (cgs)
  real(real64), parameter       :: mp_cgs = 1.672621777d-24
  !> Proton mass in kg (SI)
  real(real64), parameter       :: mp_si  = 1.672621777d-27
  !> Hydrogen mass in g (cgs)
  real(real64), parameter       :: mH_cgs = 1.6733d-24
  !> Hydrogen mass in kg (SI)
  real(real64), parameter       :: mH_si  = 1.6733d-27
  !> Electron mass in g (cgs)
  real(real64), parameter       :: me_cgs = 9.1094d-28
  !> Electron mass in kg (SI)
  real(real64), parameter       :: me_si  = 9.1094d-31
  !> Elementary charge in statcoul (cgs)
  real(real64), parameter       :: ec_cgs = 4.8032d-10
  !> Elementary charge in C (SI)
  real(real64), parameter       :: ec_si  = 1.6022d-19
  !> Boltzmann constant in erg/K (cgs)
  real(real64), parameter       :: kB_cgs = 1.3806488d-16
  !> Boltzmann constant in J/K (SI)
  real(real64), parameter       :: kB_si  = 1.3806488d-23
  !> Magnetic constant in H/m (SI)
  real(real64), parameter       :: mu0_si = 1.2566370614d-6
  !> Permittivity of free space in F/m (SI)
  real(real64), parameter       :: e0_si  = 8.8542d-12
  !> Degree of ionization (assumed fully ionized)
  real(real64), parameter       :: Z_ion = 1.0d0

contains

  !> Sets a new unit length based on the given value.
  !! @param[in] new_unit_length   New value for the unit length
  subroutine set_unit_length(new_unit_length)
    real(real64), intent(in) :: new_unit_length

    unit_length = new_unit_length

  end subroutine set_unit_length

  !> Sets a new unit numberdensity based on the given value.
  !! @param[in] new_unit_numberdensity  New value for the unit numberdensity
  subroutine set_unit_numberdensity(new_unit_numberdensity)
    real(real64), intent(in) :: new_unit_numberdensity

    unit_numberdensity = new_unit_numberdensity

  end subroutine set_unit_numberdensity

  !> Sets a new unit temperature based on the given value
  !! @param[in] new_unit_temperature  New value for the unit temperature
  subroutine set_unit_temperature(new_unit_temperature)
    real(real64), intent(in) :: new_unit_temperature

    unit_temperature = new_unit_temperature

  end subroutine set_unit_temperature

  !> Sets a new unit velocity based on the given value
  !! @param[in] new_unit_velocity   New value for the unit velocity
  subroutine set_unit_velocity(new_unit_velocity)
    real(real64), intent(in) :: new_unit_velocity

    unit_velocity = new_unit_velocity

  end subroutine set_unit_velocity

  !> Sets a new unit resistivity based on the given value. This is done in the
  !! resistivity module.
  !! @param[in] new_unit_resistivity    New value for the unit resistivity
  subroutine set_unit_resistivity(new_unit_resistivity)
    real(real64), intent(in) :: new_unit_resistivity

    unit_resistivity = new_unit_resistivity
  end subroutine set_unit_resistivity


  !> Calculates the other unit normalisations based on the pre-defined
  !! unit_length, unit_numberdensity and unit_temperature OR unit_velocity.
  !! @param[in] cgs_units   Whether to use cgs units or not.
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

    unit_dldt          = unit_luminosity / unit_temperature

    unit_conduction    = unit_density * unit_length * unit_velocity**3 / &
                         unit_temperature

    unit_dtc_drho      = unit_conduction / unit_density
    unit_dtc_dT        = unit_conduction / unit_temperature
    unit_dtc_dB2       = unit_conduction / (unit_magneticfield**2)

  end subroutine set_normalisations





end module mod_physical_constants
