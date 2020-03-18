!
! MODULE: mod_units
!
! DESCRIPTION:
!> Module containing all unit normalisations.
!
module mod_units
  use mod_global_variables, only: dp, cgs_units
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
  !> Scaling factor for gravitational acceleration
  real(dp), protected   :: unit_acceleration = 1.0d0

  !> Boolean value to check if normalisations are set
  logical, protected    :: normalisations_set = .false.

  public  :: check_if_normalisations_set
  public  :: define_nb_vel_len
  public  :: define_nb_temp_len
  public  :: define_rho_mag_len
  public  :: define_rho_temp_len
  public  :: set_unit_resistivity

  public  :: unit_length, unit_time, unit_density, unit_velocity
  public  :: unit_temperature, unit_pressure, unit_magneticfield
  public  :: unit_numberdensity, unit_luminosity, unit_dldt
  public  :: unit_conduction, unit_dtc_drho, unit_dtc_dT, unit_dtc_dB2
  public  :: unit_resistivity, unit_deta_dT, unit_acceleration

contains


  !> Checks if normalisations were set. If not, define them using default values, taken
  !! to be unit_numberdensity=1.0d9, unit_temperature=1.0d6 and unit_length=1.0d6 in cgs units.
  subroutine check_if_normalisations_set()
    if (normalisations_set) then
      return
    else
      cgs_units = .true.
      call define_nb_temp_len(new_unit_numberdensity=1.0d9, &
                              new_unit_temperature=1.0d6, &
                              new_unit_length=1.0d6)
    end if
  end subroutine check_if_normalisations_set


  !> Returns the Boltzmann constant, proton mass and magnetic permeability either in SI or cgs,
  !! based on the value of use_cgs.
  !! @param[out]  kB  Boltzmann constant
  !! @param[out]  mp  proton mass
  !! @param[out]  mu0 magnetic permeability
  subroutine get_constants(kB, mp, mu0)
    use mod_physical_constants, only: kB_si, kB_cgs, mp_si, mp_cgs, mu0_si, mu0_cgs

    real(dp), intent(out) :: kB, mp, mu0

    if (cgs_units) then
      kB = kB_cgs
      mp = mp_cgs
      mu0 = mu0_cgs
    else
      kB = kB_si
      mp = mp_si
      mu0 = mu0_si
    end if
  end subroutine get_constants


  !> Define normalisations based on a numberdensity, temperature and length unit.
  !! @param[in] new_unit_numberdensity  numberdensity
  !! @param[in] new_unit_temperature    temperature
  !! @param[in] new_unit_length         length
  subroutine define_nb_temp_len(new_unit_numberdensity, new_unit_temperature, new_unit_length)
    use mod_physical_constants, only: He_abundance

    real(dp), intent(in)  :: new_unit_numberdensity, new_unit_temperature, new_unit_length
    real(dp)  :: kB, mp, mu0

    call get_constants(kB, mp, mu0)

    unit_numberdensity = new_unit_numberdensity
    unit_temperature = new_unit_temperature
    unit_length = new_unit_length

    unit_density = (1.0d0 + 4.0d0*He_abundance) * mp * unit_numberdensity
    unit_pressure = (2.0d0 + 3.0d0*He_abundance) * unit_numberdensity * kB * unit_temperature
    unit_velocity = sqrt(unit_pressure / unit_density)
    unit_magneticfield = sqrt(mu0 * unit_pressure)

    call calculate_remaining_normalisations()
  end subroutine define_nb_temp_len


  !> Define normalisations based on a numberdensity, velocity and length unit.
  !! @param[in] new_unit_numberdensity  numberdensity
  !! @param[in] new_unit_velocity       velocity
  !! @param[in] new_unit_length         length
  subroutine define_nb_vel_len(new_unit_numberdensity, new_unit_velocity, new_unit_length)
    use mod_physical_constants, only: He_abundance

    real(dp), intent(in)  :: new_unit_numberdensity, new_unit_velocity, new_unit_length
    real(dp)  :: kB, mp, mu0

    call get_constants(kB, mp, mu0)

    unit_numberdensity = new_unit_numberdensity
    unit_velocity = new_unit_velocity
    unit_length = new_unit_length

    unit_density = (1.0d0 + 4.0d0*He_abundance) * mp * unit_numberdensity
    unit_pressure = unit_density * unit_velocity**2
    unit_temperature = unit_pressure / ((2.0d0 + 3.0d0*He_abundance) * unit_numberdensity * kB)
    unit_magneticfield = sqrt(mu0 * unit_pressure)

    call calculate_remaining_normalisations()
  end subroutine define_nb_vel_len


  !> Define normalisations based on a density, temperature and length unit.
  !! @param[in] new_unit_density        density
  !! @param[in] new_unit_temperature    temperature
  !! @param[in] new_unit_length         length
  subroutine define_rho_temp_len(new_unit_density, new_unit_temperature, new_unit_length)
    use mod_physical_constants, only: He_abundance

    real(dp), intent(in)  :: new_unit_density, new_unit_temperature, new_unit_length
    real(dp)  :: kB, mp, mu0

    call get_constants(kB, mp, mu0)

    unit_density = new_unit_density
    unit_temperature = new_unit_temperature
    unit_length = new_unit_length

    unit_numberdensity = unit_density / ((1.0d0 + 4.0d0*He_abundance) * mp)
    unit_pressure = (2.0d0 + 3.0d0*He_abundance) * unit_numberdensity * kB * unit_temperature
    unit_velocity = sqrt(unit_pressure / unit_density)
    unit_magneticfield = sqrt(mu0 * unit_pressure)

    call calculate_remaining_normalisations()
  end subroutine define_rho_temp_len


  !> Define normalisations based on a density, magnetic field and length unit.
  !! @param[in] new_unit_density        density
  !! @param[in] new_unit_magneticfield  magnetic field
  !! @param[in] new_unit_length         length
  subroutine define_rho_mag_len(new_unit_density, new_unit_magneticfield, new_unit_length)
    use mod_physical_constants, only: He_abundance

    real(dp), intent(in)  :: new_unit_density, new_unit_magneticfield, new_unit_length
    real(dp)  :: kB, mp, mu0

    call get_constants(kB, mp, mu0)

    unit_density = new_unit_density
    unit_magneticfield = new_unit_magneticfield
    unit_length = new_unit_length

    unit_numberdensity = unit_density / ((1.0d0 + 4.0d0*He_abundance) * mp)
    unit_pressure = unit_magneticfield**2 / mu0
    unit_temperature = unit_pressure / ((2.0d0 + 3.0d0*He_abundance) * unit_numberdensity * kB)
    unit_velocity = sqrt(unit_pressure / unit_density)

    call calculate_remaining_normalisations()
  end subroutine define_rho_mag_len


  !> Calculates the remaining normalisations which are based on density, magneticfield, length etc.
  !! These have been set at this point.
  subroutine calculate_remaining_normalisations()
    use mod_physical_constants, only: He_abundance

    unit_time          = unit_length / unit_velocity

    unit_luminosity    = unit_pressure / ((1.0d0 + 2.0d0*He_abundance) * (unit_numberdensity**2 * unit_time))
    unit_dldt          = unit_luminosity / unit_temperature

    unit_conduction    = unit_density * unit_length * unit_velocity**3 / unit_temperature
    unit_dtc_drho      = unit_conduction / unit_density
    unit_dtc_dT        = unit_conduction / unit_temperature
    unit_dtc_dB2       = unit_conduction / (unit_magneticfield**2)

    unit_acceleration  = unit_length / unit_time**2

    normalisations_set = .true.
  end subroutine calculate_remaining_normalisations


  !> Sets the unit resistivity based on the given value. This is done automatically in the
  !! resistivity module, and is set as such that the normalised resistivity at 1 MK equals 0.1.
  !! @param[in] new_unit_resistivity  resistivity
  subroutine set_unit_resistivity(new_unit_resistivity)
    real(dp), intent(in)  :: new_unit_resistivity

    unit_resistivity = new_unit_resistivity
    unit_deta_dT = unit_resistivity / unit_temperature
  end subroutine set_unit_resistivity

end module mod_units
