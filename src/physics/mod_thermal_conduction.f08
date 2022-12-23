! =============================================================================
!> This module is responsible for calculating and setting the
!! thermal conduction values based on the equilibrium configuration.
module mod_thermal_conduction
  use mod_global_variables, only: dp
  use mod_physical_constants, only: dpi, coulomb_log
  use mod_logging, only: log_message, str
  use mod_settings, only: settings_t
  implicit none

  private

  real(dp)  :: pf_kappa_para
  real(dp)  :: pf_kappa_perp

  public :: set_conduction_values

contains


  !> This routines sets all thermal conduction values in <tt>kappa_field</tt>,
  !! and calls all other relevant subroutines defined in this module.
  !! @note    The derivatives of the perpendicular thermal conduction terms
  !!          are all zero if that value is taken to be fixed.
  subroutine set_conduction_values(settings, rho_field, T_field, B_field, kappa_field)
    use mod_types, only: density_type, temperature_type, bfield_type, conduction_type

    type(settings_t), intent(in) :: settings
    !> the type containing the density attributes
    type(density_type), intent(in)        :: rho_field
    !> the type containing the temperature attributes
    type(temperature_type), intent(in)    :: T_field
    !> the type containing the magnetic field attributes
    type(bfield_type), intent(in)         :: B_field
    !> the type containing the thermal conduction attributes
    type(conduction_type), intent(inout)  :: kappa_field

    ! set conduction prefactor in cgs units
    pf_kappa_para = 1.8d-5
    pf_kappa_perp = 8.2d-10

    if (settings%physics%conduction%has_parallel_conduction()) then
      call set_kappa_para(settings, T_field%T0, kappa_field%kappa_para)
      call set_kappa_para_derivatives(settings, T_field%T0, kappa_field%d_kappa_para_dT)
    end if
    if (settings%physics%conduction%has_perpendicular_conduction()) then
      call set_kappa_perp( &
        settings, T_field%T0, rho_field%rho0, B_field%B0, kappa_field%kappa_perp &
      )
      ! set rho, T, B derivatives
      call set_kappa_perp_derivatives( &
        settings, &
        T_field % T0, &
        rho_field % rho0, &
        B_field % B0, &
        kappa_field % d_kappa_perp_drho, &
        kappa_field % d_kappa_perp_dB2, &
        kappa_field % d_kappa_perp_dT &
      )
      ! set derivatives with respect to r
      call set_kappa_perp_radial_derivative( &
        settings, rho_field, T_field, B_field, kappa_field &
      )
    end if

    ! set conduction prefactor and its radial derivative
    call set_conduction_prefactor(T_field, B_field, kappa_field)
  end subroutine set_conduction_values


  !> Calculates the parallel thermal conduction.
  !! Returns either the full parallel thermal conduction based on
  !! the equilibrium parameters, or a fixed value if specified
  !! in the global variables module.
  !! @note    The temperature should be normalised when calling this routine,
  !!          the parallel conduction is normalised on exit.
  subroutine set_kappa_para(settings, T0_eq, tc_para)
    type(settings_t), intent(in) :: settings
    !> equilibrium temperature
    real(dp), intent(in)  :: T0_eq(:)
    !> parallel thermal conduction values, normalised on exit
    real(dp), intent(out) :: tc_para(size(T0_eq))
    real(dp)              :: T0_dimfull(size(T0_eq))

    if (settings%physics%conduction%has_fixed_tc_para()) then
      tc_para = settings%physics%conduction%get_fixed_tc_para()
      return
    end if

    T0_dimfull = T0_eq * settings%units%get_unit_temperature()
    tc_para = pf_kappa_para * T0_dimfull**2.5d0 / coulomb_log
    tc_para = tc_para / settings%units%get_unit_conduction()
  end subroutine set_kappa_para


  !> Calculates the perpendicular thermal conduction.
  !! Returns either the full perpendicular thermal conduction based on
  !! the equilibrium parameters, or a fixed value if specified
  !! in the global variables module.
  !! @note    All variables should be normalised when calling this routine,
  !!          the perpendicular conduction is normalised on exit.
  subroutine set_kappa_perp(settings, T0_eq, rho0_eq, B0_eq, tc_perp)
    type(settings_t), intent(in) :: settings
    !> equilibrium temperature
    real(dp), intent(in)  :: T0_eq(:)
    !> equilibrium density
    real(dp), intent(in)  :: rho0_eq(:)
    !> equilibrium magnetic field
    real(dp), intent(in)  :: B0_eq(:)
    !> perpendicular thermal conduction values, normalised on exit
    real(dp), intent(out) :: tc_perp(size(T0_eq))
    real(dp) :: T0_dimfull(size(T0_eq)), B0_dimfull(size(T0_eq))
    real(dp) :: nH_dimfull(size(T0_eq))

    if (settings%physics%conduction%has_fixed_tc_perp()) then
      tc_perp = settings%physics%conduction%get_fixed_tc_perp()
      return
    end if

    nH_dimfull = rho0_eq * settings%units%get_unit_numberdensity()
    T0_dimfull = T0_eq * settings%units%get_unit_temperature()
    B0_dimfull = B0_eq * settings%units%get_unit_magneticfield()

    tc_perp = pf_kappa_para * pf_kappa_perp * coulomb_log * nH_dimfull**2 &
      / (B0_dimfull**2 * sqrt(T0_dimfull))
    tc_perp = tc_perp / settings%units%get_unit_conduction()
  end subroutine set_kappa_perp


  !> Calculates the temperature derivative of the parallel thermal conduction component.
  !! @note    All variables should be normalised when calling this routine,
  !!          values are normalised on exit.
  subroutine set_kappa_para_derivatives(settings, T0_eq, d_tcpara_dT)
    type(settings_t), intent(in) :: settings
    !> equilibrium temperature
    real(dp), intent(in)  :: T0_eq(:)
    !> temperature derivative of parallel thermal conduction, normalised on exit
    real(dp), intent(out) :: d_tcpara_dT(size(T0_eq))
    real(dp) :: T0_dimfull(size(T0_eq))
    real(dp) :: unit_temperature, unit_dtc_dT

    if (settings%physics%conduction%has_fixed_tc_para()) return

    call log_message("setting kappa_para derivatives", level="debug")
    unit_temperature = settings%units%get_unit_temperature()
    unit_dtc_dT = settings%units%get_unit_conduction() / unit_temperature
    T0_dimfull = T0_eq * unit_temperature
    d_tcpara_dT = pf_kappa_para * 2.5d0 * T0_dimfull**1.5d0 / coulomb_log
    d_tcpara_dT = d_tcpara_dT / unit_dtc_dT
  end subroutine set_kappa_para_derivatives


  !> Calculates the thermal conduction derivatives.
  !! Returns the derivative of the perpendicular thermal conduction with respect to
  !! density, magnetic field squared and temperature.
  !! @note    All variables should be normalised when calling this routine.
  !!          The thermal conduction derivatives are normalised on exit.
  subroutine set_kappa_perp_derivatives( &
    settings, T0_eq, rho0_eq, B0_eq, d_tc_drho, d_tc_dB2, d_tc_dT &
  )
    type(settings_t), intent(in) :: settings
    !> equilibrium temperature
    real(dp), intent(in)  :: T0_eq(:)
    !> equilibrium density
    real(dp), intent(in)  :: rho0_eq(:)
    !> equilibrium magnetic field
    real(dp), intent(in)  :: B0_eq(:)
    !> derivative of perpendicular conduction with respect to density
    real(dp), intent(out) :: d_tc_drho(size(T0_eq))
    !> derivative of perpendicular conduction with respect to magnetic field squared
    real(dp), intent(out) :: d_tc_dB2(size(T0_eq))
    !> derivative of perpendicular conduction with respect to temperature
    real(dp), intent(out) :: d_tc_dT(size(T0_eq))
    real(dp)  :: T0_dimfull(size(T0_eq))
    real(dp)  :: rho0_dimfull(size(T0_eq))
    real(dp)  :: B0_dimfull(size(T0_eq))
    real(dp)  :: nH_dimfull(size(T0_eq))
    real(dp) :: unit_conduction

    if (settings%physics%conduction%has_fixed_tc_perp()) return

    call log_message("setting kappa_perp rho, T, B derivatives", level="debug")
    unit_conduction = settings%units%get_unit_conduction()
    ! re-dimensionalise variables, in normalised units nH = rho
    nH_dimfull = rho0_eq * settings%units%get_unit_numberdensity()
    T0_dimfull = T0_eq * settings%units%get_unit_temperature()
    rho0_dimfull = rho0_eq * settings%units%get_unit_density()
    B0_dimfull = B0_eq * settings%units%get_unit_magneticfield()

    ! density derivative
    d_tc_drho = 2.0d0 * pf_kappa_para * pf_kappa_perp * coulomb_log * nH_dimfull &
      / (B0_dimfull**2 * sqrt(T0_dimfull))
    ! magnetic field derivative
    d_tc_dB2 = -pf_kappa_para * pf_kappa_perp * coulomb_log * nH_dimfull**2 &
      / (B0_dimfull**4 * sqrt(T0_dimfull))
    ! temperature derivative
    d_tc_dT = -0.5d0 * pf_kappa_para * pf_kappa_perp * coulomb_log * nH_dimfull**2 &
      / (B0_dimfull**2 * T0_dimfull**(3.0d0/2.0d0))
    d_tc_drho = d_tc_drho / (unit_conduction / settings%units%get_unit_density())
    d_tc_dB2 = d_tc_dB2 / (unit_conduction / settings%units%get_unit_magneticfield()**2)
    d_tc_dT = d_tc_dT / (unit_conduction / settings%units%get_unit_temperature())
  end subroutine set_kappa_perp_derivatives


  !> Sets the radial derivative of the perpendicular thermal conduction coefficient.
  !! This is defined as
  !! $$
  !! kappa_{\perp, 0}' = \frac{d\kappa_{\perp, 0}}{d\rho} \rho_0'
  !!  + \frac{d\kappa_{\perp, 0}}{dT} T_0'
  !!  + \frac{d\kappa_{\perp, 0}}{d(B^2)}2B_0 B_0'
  !! $$
  subroutine set_kappa_perp_radial_derivative( &
    settings, rho_field, T_field, B_field, kappa_field &
  )
    use mod_types, only: density_type, temperature_type, bfield_type, conduction_type

    type(settings_t), intent(in) :: settings
    !> the type containing the density attributes
    type(density_type), intent(in) :: rho_field
    !> the type containing the temperature attributes
    type(temperature_type), intent(in)  :: T_field
    !> the type containing the magnetic field attributes
    type(bfield_type), intent(in) :: B_field
    !> the type containing the thermal conduction attributes
    type(conduction_type), intent(inout)  :: kappa_field

    real(dp)  :: dB0(size(B_field % B0))

    if (settings%physics%conduction%has_fixed_tc_perp()) return

    call log_message("setting kappa_perp radial derivative", level="debug")
    ! magnetic field derivative
    dB0 = ( &
      B_field % B02 * B_field % d_B02_dr + B_field % B03 * B_field % d_B03_dr &
    ) / B_field % B0
    ! coordinate derivative of kappa_perp
    kappa_field % d_kappa_perp_dr = ( &
      kappa_field % d_kappa_perp_drho * rho_field % d_rho0_dr &
      + kappa_field % d_kappa_perp_dT * T_field % d_T0_dr &
      + kappa_field % d_kappa_perp_dB2 * 2.0d0 * B_field % B0 * dB0 &
    )
  end subroutine set_kappa_perp_radial_derivative


  !> Sets the thermal conduction prefactor, given by
  !! $$ \frac{\kappa_{\parallel, 0} - \kappa_{\perp, 0}}{B_0^2} $$.
  !! The radial derivative of the prefactor is also set, given by
  !! $$ \frac{1}{B_0^3}\left[
  !!  -2\left(\kappa_{\parallel, 0} - \kappa_{\perp, 0}\right)B_0'
  !!  +\left(\kappa_{\parallel, 0}' - \kappa_{\perp, 0}'\right)B_0
  !! \right]
  !! $$
  subroutine set_conduction_prefactor(T_field, B_field, kappa_field)
    use mod_types, only: temperature_type, bfield_type, conduction_type

    !> the type containing the temperature attributes
    type(temperature_type), intent(in)  :: T_field
    !> the type containing the magnetic field attributes
    type(bfield_type), intent(in) :: B_field
    !> the type containing the thermal conduction attributes
    type(conduction_type), intent(inout)  :: kappa_field

    real(dp)  :: d_kappa_para_dr(size(T_field % T0))
    real(dp)  :: dB0(size(B_field % B0))

    call log_message("setting conduction prefactor", level="debug")
    ! calculate and set conduction prefactor
    kappa_field % prefactor = ( &
      (kappa_field % kappa_para - kappa_field % kappa_perp) / (B_field % B0**2) &
    )

    ! radial derivative of parallel conduction component
    d_kappa_para_dr = kappa_field % d_kappa_para_dT * T_field % d_T0_dr
    ! radial derivative of magnetic field
    dB0 = ( &
      B_field % B02 * B_field % d_B02_dr + B_field % B03 * B_field % d_B03_dr &
    ) / B_field % B0

    call log_message("setting conduction prefactor radial derivative", level="debug")
    kappa_field % d_prefactor_dr = ( &
      (d_kappa_para_dr - kappa_field % d_kappa_perp_dr) * B_field % B0 &
      - 2.0d0 * (kappa_field % kappa_para - kappa_field % kappa_perp) * dB0 &
    ) / (B_field % B0**3)
  end subroutine set_conduction_prefactor

end module mod_thermal_conduction
