! =============================================================================
!> This module is responsible for calculating and setting the
!! thermal conduction values based on the equilibrium configuration.
module mod_thermal_conduction
  use mod_global_variables, only: dp
  use mod_physical_constants, only: dpi, coulomb_log, tc_pf_kappa_para, tc_pf_kappa_perp
  use mod_logging, only: logger, str
  use mod_settings, only: settings_t
  use mod_check_values, only: is_zero
  use mod_background, only: background_t
  use mod_function_utils, only: from_function
  use mod_grid, only: grid_gauss
  use mod_types, only: conduction_type
  implicit none

  private

  public :: set_conduction_values

contains


  !> This routines sets all thermal conduction values in <tt>kappa_field</tt>,
  !! and calls all other relevant subroutines defined in this module.
  !! @note    The derivatives of the perpendicular thermal conduction terms
  !!          are all zero if that value is taken to be fixed.
  subroutine set_conduction_values(settings, background, kappa_field)
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    !> the type containing the thermal conduction attributes
    type(conduction_type), intent(inout)  :: kappa_field

    if (.not. settings%has_bfield()) call logger%info( &
      "no B-field detected, using isotropic thermal conduction" &
    )

    call logger%debug("setting thermal conduction values")
    if (settings%physics%conduction%has_parallel_conduction()) then
      call set_kappa_para( &
        settings, from_function(background%temperature%T0, grid_gauss), &
        kappa_field%kappa_para &
      )
      call set_kappa_para_derivatives( &
        settings, &
        from_function(background%temperature%T0, grid_gauss), &
        kappa_field%d_kappa_para_dT &
      )
    end if
    if (settings%physics%conduction%has_perpendicular_conduction()) then
      call set_kappa_perp( &
        settings=settings, &
        T0_eq=from_function(background%temperature%T0, grid_gauss), &
        rho0_eq=from_function(background%density%rho0, grid_gauss), &
        B0_eq=background%magnetic%get_B0(grid_gauss), &
        tc_perp=kappa_field%kappa_perp &
      )
      ! set rho, T, B derivatives
      call set_kappa_perp_derivatives( &
        settings=settings, &
        T0_eq=from_function(background%temperature%T0, grid_gauss), &
        rho0_eq=from_function(background%density%rho0, grid_gauss), &
        B0_eq=background%magnetic%get_B0(grid_gauss), &
        d_tc_drho=kappa_field % d_kappa_perp_drho, &
        d_tc_dB2=kappa_field % d_kappa_perp_dB2, &
        d_tc_dT=kappa_field % d_kappa_perp_dT &
      )
      ! set derivatives with respect to r
      call set_kappa_perp_radial_derivative(settings, background, kappa_field)
    end if

    ! set conduction prefactor and its radial derivative
    call set_conduction_prefactor(settings, background, kappa_field)
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
      call logger%debug("  using fixed parallel thermal conduction")
      tc_para = settings%physics%conduction%get_fixed_tc_para()
      return
    end if

    call logger%debug("  using T-dependent parallel thermal conduction")
    T0_dimfull = T0_eq * settings%units%get_unit_temperature()
    tc_para = tc_pf_kappa_para * T0_dimfull**2.5d0 / coulomb_log
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
      call logger%debug("  using fixed perpendicular thermal conduction")
      tc_perp = settings%physics%conduction%get_fixed_tc_perp()
      return
    end if
    if (.not. settings%has_bfield()) then
      call logger%debug("  no B-field, kappa_para = kappa_perp")
      call set_kappa_para(settings, T0_eq, tc_perp)
      return
    end if

    call logger%debug("  using (rho, T, B)-dependent perpendicular thermal conduction")
    nH_dimfull = rho0_eq * settings%units%get_unit_numberdensity()
    T0_dimfull = T0_eq * settings%units%get_unit_temperature()
    B0_dimfull = B0_eq * settings%units%get_unit_magneticfield()

    tc_perp = tc_pf_kappa_para * tc_pf_kappa_perp * coulomb_log * nH_dimfull**2 &
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

    unit_temperature = settings%units%get_unit_temperature()
    unit_dtc_dT = settings%units%get_unit_conduction() / unit_temperature
    T0_dimfull = T0_eq * unit_temperature
    d_tcpara_dT = tc_pf_kappa_para * 2.5d0 * T0_dimfull**1.5d0 / coulomb_log
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

    call logger%debug("  setting kappa_perp derivatives")
    if (.not. settings%has_bfield()) then
      call logger%debug("  no Bfield, dkappa_perp/dT = dkappa_para/dT, others zero")
      call set_kappa_para_derivatives(settings, T0_eq, d_tc_dT)
      d_tc_drho = 0.0_dp
      d_tc_dB2 = 0.0_dp
      return
    end if

    unit_conduction = settings%units%get_unit_conduction()
    ! re-dimensionalise variables, in normalised units nH = rho
    nH_dimfull = rho0_eq * settings%units%get_unit_numberdensity()
    T0_dimfull = T0_eq * settings%units%get_unit_temperature()
    rho0_dimfull = rho0_eq * settings%units%get_unit_density()
    B0_dimfull = B0_eq * settings%units%get_unit_magneticfield()

    ! density derivative
    d_tc_drho = 2.0d0 * tc_pf_kappa_para * tc_pf_kappa_perp * coulomb_log * nH_dimfull &
      / (B0_dimfull**2 * sqrt(T0_dimfull))
    ! magnetic field derivative
    d_tc_dB2 = -tc_pf_kappa_para * tc_pf_kappa_perp * coulomb_log * nH_dimfull**2 &
      / (B0_dimfull**4 * sqrt(T0_dimfull))
    ! temperature derivative
    d_tc_dT = ( &
      -0.5d0 * tc_pf_kappa_para * tc_pf_kappa_perp * coulomb_log * nH_dimfull**2 &
    ) / (B0_dimfull**2 * T0_dimfull**(3.0d0/2.0d0))
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
  subroutine set_kappa_perp_radial_derivative(settings, background, kappa_field)
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    !> the type containing the thermal conduction attributes
    type(conduction_type), intent(inout)  :: kappa_field

    real(dp) :: drho0(settings%grid%get_gauss_gridpts())
    real(dp) :: dT0(settings%grid%get_gauss_gridpts())
    real(dp) :: B0(settings%grid%get_gauss_gridpts())
    real(dp) :: dB0(settings%grid%get_gauss_gridpts())

    if (settings%physics%conduction%has_fixed_tc_perp()) return

    call logger%debug("  setting kappa_perp radial derivative")
    if (.not. settings%has_bfield()) then
      call logger%debug("  no B-field, derivative is dkappa_para/dT * T0'")
      kappa_field % d_kappa_perp_dr = ( &
        kappa_field % d_kappa_perp_dT &
        * from_function(background%temperature%dT0, grid_gauss) &
      )
      return
    end if
    drho0 = from_function(background%density%drho0, grid_gauss)
    dT0 = from_function(background%temperature%dT0, grid_gauss)
    B0 = background%magnetic%get_B0(grid_gauss)
    dB0 = background%magnetic%get_dB0(grid_gauss)
    ! coordinate derivative of kappa_perp
    kappa_field % d_kappa_perp_dr = ( &
      kappa_field % d_kappa_perp_drho * drho0 &
      + kappa_field % d_kappa_perp_dT * dT0 &
      + kappa_field % d_kappa_perp_dB2 * 2.0d0 * B0 * dB0 &
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
  subroutine set_conduction_prefactor(settings, background, kappa_field)
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    !> the type containing the thermal conduction attributes
    type(conduction_type), intent(inout)  :: kappa_field

    real(dp) :: d_kappa_para_dr(settings%grid%get_gauss_gridpts())
    real(dp) :: drho0(settings%grid%get_gauss_gridpts())
    real(dp) :: dT0(settings%grid%get_gauss_gridpts())
    real(dp) :: B0(settings%grid%get_gauss_gridpts())
    real(dp) :: dB0(settings%grid%get_gauss_gridpts())

    call logger%debug("  setting conduction prefactor")
    if (.not. settings%has_bfield()) then
      call logger%debug("  no B-field, (kappa_para - kappa_perp) is zero")
      kappa_field%prefactor = 0.0_dp
      kappa_field%d_prefactor_dr = 0.0_dp
      return
    end if

    drho0 = from_function(background%density%drho0, grid_gauss)
    dT0 = from_function(background%temperature%dT0, grid_gauss)
    B0 = background%magnetic%get_B0(grid_gauss)
    dB0 = background%magnetic%get_dB0(grid_gauss)

    ! calculate and set conduction prefactor
    kappa_field % prefactor = ( &
      (kappa_field % kappa_para - kappa_field % kappa_perp) / (B0**2) &
    )
    ! radial derivative of parallel conduction component
    d_kappa_para_dr = kappa_field % d_kappa_para_dT * dT0

    call logger%debug("  setting conduction prefactor radial derivative")
    kappa_field % d_prefactor_dr = ( &
      (d_kappa_para_dr - kappa_field % d_kappa_perp_dr) * B0 &
      - 2.0d0 * (kappa_field % kappa_para - kappa_field % kappa_perp) * dB0 &
    ) / B0**3
  end subroutine set_conduction_prefactor

end module mod_thermal_conduction
