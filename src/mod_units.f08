! =============================================================================
!> This module handles checking and setting all unit normalisations.
!! When setting the normalisations, one should always specify a unit length
!! and a unit magneticfield. The normalisations use a reference value of 2 for plasma beta
!! as is common practice, such that we can write \(\beta/2 = \mu * p / B^2 = 1 \).
!! The pressure unit is then calculated as
!! $$ p_{unit} = \frac{B_{unit}^2}{\mu} $$
!! Then, if a unit density is specified this fixes the unit temperature and vice-versa,
!! through use of the ideal gas law
!! $$ p_{unit} = \mathcal{R}_{specific}T_{unit}\rho_{unit} $$
!! The other normalisations are then fixed, and are set through
!! $$ n_{unit} = \frac{\rho_{unit}}{m_p}, ~~~~~~~
!!    v_{unit} = \frac{B_{unit}}{\sqrt{\mu\rho_{unit}}},  $$
!! $$ t_{unit} = \frac{L_{unit}}{v_{unit}},  ~~~~~~~~~
!!    \mathscr{L}_{unit} = \frac{p_{unit}}{t_{unit}n_{unit}^2}, $$
!! $$ \kappa_{unit} = \frac{\rho_{unit}L_{unit}v_{unit}^3}{T_{unit}} $$
!! Where the latter two denote the luminosity and thermal conduction normalisations.
!! @note    The normalisations are only used when conduction, resistivity or radiative cooling
!!          is included. However, we always set base units for good practice, and take reference
!!          values of 1 MK, 10 Gauss and <tt>1e9</tt> cm. @endnote
!! @note    As one can always specify a unit current, the unit resistivity is defined in such
!!          a way that the normalised resistivity equals 0.1 at 1 MK. @endnote
module mod_units
  use mod_global_variables, only: dp, cgs_units
  use mod_logging, only: log_message, str
  implicit none

  private

  !> length normalisation
  real(dp), protected   :: unit_length = 1.0d0
  !> time normalisation
  real(dp), protected   :: unit_time = 1.0d0
  !> density normalisation
  real(dp), protected   :: unit_density = 1.0d0
  !> velocity normalisation
  real(dp), protected   :: unit_velocity = 0.0d0
  !> temperature normalisation
  real(dp), protected   :: unit_temperature = 1.0d0
  !> pressure normalisation
  real(dp), protected   :: unit_pressure = 1.0d0
  !> magnetic field normalisation
  real(dp), protected   :: unit_magneticfield = 1.0d0
  !> numberdensity normalisation
  real(dp), protected   :: unit_numberdensity = 1.0d0
  !> cooling function normalisation
  real(dp), protected   :: unit_lambdaT = 1.0d0
  !> normalisation for derivative of cooling function with respect to temperature
  real(dp), protected   :: unit_dlambdaT_dT = 1.0d0
  !> thermal conduction normalisation
  real(dp), protected   :: unit_conduction = 1.0d0
  !> normalisation for d(kappa_perp)/d(rho)
  real(dp), protected   :: unit_dtc_drho = 1.0d0
  !> normalisation for d(kappa_perp)/d(T)
  real(dp), protected   :: unit_dtc_dT = 1.0d0
  !> normalisation for d(kappa_perp)/d(B^2)
  real(dp), protected   :: unit_dtc_dB2 = 1.0d0
  !> resistivity normalisation
  real(dp), protected   :: unit_resistivity = 1.0d0
  !> normalisation for d(eta)/d(T)
  real(dp), protected   :: unit_deta_dT = 1.0d0
  !> mass normalisation
  real(dp), protected   :: unit_mass = 1.0d0
  !> mean molecular weight
  real(dp), protected   :: mean_molecular_weight = 1.0d0


  !> boolean to check if normalisations are set
  logical, protected    :: normalisations_are_set = .false.

  public  :: check_if_normalisations_set
  public  :: normalisations_are_set
  public  :: set_normalisations
  public  :: set_unit_resistivity

  public  :: unit_length, unit_time, unit_density, unit_velocity
  public  :: unit_temperature, unit_pressure, unit_magneticfield
  public  :: unit_numberdensity, unit_lambdaT, unit_dlambdaT_dT
  public  :: unit_conduction, unit_dtc_drho, unit_dtc_dT, unit_dtc_dB2
  public  :: unit_resistivity, unit_deta_dT
  public  :: mean_molecular_weight

contains


  !> Checks if normalisations are set.
  !! If normalisations are not set, set them to default values.
  !! These are 1 MK as unit temperature, 10 Gauss as unit magnetic field and
  !! <tt>1e9</tt> cm as a unit length. If normalisations are already set through
  !! the equilibrium submodule or parfiles nothing is done.
  subroutine check_if_normalisations_set()
    if (normalisations_are_set) then
      call log_message("normalisations are already set", level="debug")
      return
    else
      cgs_units = .true.
      call set_normalisations( &
        new_unit_temperature=1.0d6, &
        new_unit_magneticfield=10.0d0, &
        new_unit_length=1.0d9 &
      )
    end if
  end subroutine check_if_normalisations_set


  !> Returns the Boltzmann constant, proton mass, magnetic constant
  !! and gas constant in either cgs (default) or SI units.
  subroutine get_constants(kB, mp, mu0, Rgas)
    use mod_physical_constants, only: kB_si, kB_cgs, mp_si, mp_cgs, mu0_si, mu0_cgs, R_si, R_cgs

    !> Boltzmann constant
    real(dp), intent(out) :: kB
    !> proton mass
    real(dp), intent(out) :: mp
    !> magnetic constant
    real(dp), intent(out) :: mu0
    !> gas constant
    real(dp), intent(out) :: Rgas

    if (cgs_units) then
      call log_message("getting constants in cgs units", level="debug")
      kB = kB_cgs
      mp = mp_cgs
      mu0 = mu0_cgs
      Rgas = R_cgs
    else
      call log_message("getting constants in SI units", level="debug")
      kB = kB_si
      mp = mp_si
      mu0 = mu0_si
      Rgas = R_si
    end if
  end subroutine get_constants


  !> Defines unit normalisations based on a magnetic field unit, length unit,
  !! and a density OR temperature unit. Calling this routine automatically
  !! sets <tt>normalisations_are_set</tt> to <tt>True</tt>.
  !! An optional mean molecular weight can be passed, which defaults to 1/2
  !! corresponding to an electron-proton plasma.
  !! @note  It is good practice to specify the keyword arguments when calling this
  !!        routine, to make sure units are not switched. @endnote
  !! @warning Throws an error if:
  !!
  !! - the unit density and unit temperature are both specified.
  !! - neither unit density or unit temperature is specified.  @endwarning
  subroutine set_normalisations( &
    new_unit_density, &
    new_unit_temperature, &
    new_unit_magneticfield, &
    new_unit_length, &
    new_mean_molecular_weight &
  )
    !> new value for the unit density
    real(dp), intent(in), optional  :: new_unit_density
    !> new value for the unit temperature
    real(dp), intent(in), optional  :: new_unit_temperature
    !> new value for the unit magnetic field
    real(dp), intent(in)            :: new_unit_magneticfield
    !> new value for the unit length
    real(dp), intent(in)            :: new_unit_length
    !> new value for the mean molecular weigth
    real(dp), intent(in), optional  :: new_mean_molecular_weight
    real(dp)  :: kB, mp, mu0, Rgas

    call get_constants(kB, mp, mu0, Rgas)

    if (present(new_unit_density) .and. present(new_unit_temperature)) then
      call log_message( &
        "unit density and unit temperature can not both be set.", level="error" &
      )
    end if

    if (present(new_mean_molecular_weight)) then
      mean_molecular_weight = new_mean_molecular_weight
      call log_message( &
        "mean molecular weight set to " // str(mean_molecular_weight), level="info" &
      )
    end if

    ! TODO (niels):
    ! remove this error once SI is properly checked. I suspect there is an
    ! inconsistency in the radiative cooling module: all tables are given in
    ! cgs units, so I think if SI units are specified we first have to scale
    ! the tables to SI and THEN normalise using SI normalisations.
    if (.not. cgs_units) then
      call log_message( &
        "possible inconsistency in SI units, use cgs for now!", level="error" &
      )
    end if

    unit_magneticfield = new_unit_magneticfield
    unit_length = new_unit_length
    unit_pressure = unit_magneticfield**2 / mu0

    if (present(new_unit_density)) then
      unit_density = new_unit_density
      unit_temperature = ( &
        mean_molecular_weight * unit_pressure * mp / (kB * unit_density) &
      )
    else if (present(new_unit_temperature)) then
      unit_temperature = new_unit_temperature
      unit_density = ( &
        mean_molecular_weight * unit_pressure * mp / (kB * unit_temperature) &
      )
    else
      call log_message("no unit density or unit temperature specified.", level="error")
    end if

    unit_mass = unit_density * unit_length**3
    unit_numberdensity = unit_density / mp
    unit_velocity = unit_magneticfield / sqrt(mu0 * unit_density)
    unit_time = unit_length / unit_velocity
    unit_lambdaT = unit_pressure / (unit_time * unit_numberdensity**2)
    unit_dlambdaT_dT = unit_lambdaT / unit_temperature

    unit_conduction    = unit_density * unit_length * unit_velocity**3 / unit_temperature
    unit_dtc_drho      = unit_conduction / unit_density
    unit_dtc_dT        = unit_conduction / unit_temperature
    unit_dtc_dB2       = unit_conduction / (unit_magneticfield**2)

    normalisations_are_set = .true.
  end subroutine set_normalisations


  !> Sets the unit resistivity.
  !! This routine is called by the resistivity module and is used to set the
  !! resistivity unit in such a way that eta equals 0.1 when the temperature is 1 MK.
  subroutine set_unit_resistivity(new_unit_resistivity)
    !> new value for the unit resistivity
    real(dp), intent(in)  :: new_unit_resistivity

    unit_resistivity = new_unit_resistivity
    unit_deta_dT = unit_resistivity / unit_temperature
  end subroutine set_unit_resistivity

end module mod_units
