! =============================================================================
!> This module is responsible for calculating and setting the
!! thermal conduction values based on the equilibrium configuration.
module mod_thermal_conduction
  use mod_global_variables, only: dp, gauss_gridpts, cgs_units
  use mod_physical_constants, only: dpi, coulomb_log, mp_cgs, mp_si
  use mod_units, only: unit_temperature, unit_conduction, unit_magneticfield, &
                       unit_numberdensity, unit_density
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
  subroutine set_conduction_values(rho_field, T_field, B_field, kappa_field)
    use mod_types, only: density_type, temperature_type, bfield_type, conduction_type
    use mod_global_variables, only: use_fixed_tc_para, use_fixed_tc_perp

    !> the type containing the density attributes
    type(density_type), intent(in)        :: rho_field
    !> the type containing the temperature attributes
    type(temperature_type), intent(in)    :: T_field
    !> the type containing the magnetic field attributes
    type(bfield_type), intent(in)         :: B_field
    !> the type containing the thermal conduction attributes
    type(conduction_type), intent(inout)  :: kappa_field

    ! set conduction prefactor depending on unit system
    if (cgs_units) then
      pf_kappa_para = 1.8d-5
      pf_kappa_perp = 8.2d-10
    else
      pf_kappa_para = 1.8d-10
      pf_kappa_perp = 8.2d-33
    end if

    call get_kappa_para(T_field % T0, kappa_field % kappa_para)
    call get_kappa_perp( &
      T_field % T0, rho_field % rho0, B_field % B0, kappa_field % kappa_perp &
    )

    if (.not. use_fixed_tc_para) then
      call get_kappa_para_derivatives(T_field % T0, kappa_field % d_kappa_para_dT)
    end if

    if (.not. use_fixed_tc_perp) then
      call get_kappa_perp_derivatives( &
        T_field % T0, &
        rho_field % rho0, &
        B_field % B0, &
        kappa_field % d_kappa_perp_drho, &
        kappa_field % d_kappa_perp_dB2, &
        kappa_field % d_kappa_perp_dT &
      )
    end if
  end subroutine set_conduction_values


  !> Calculates the parallel thermal conduction.
  !! Returns either the full parallel thermal conduction based on
  !! the equilibrium parameters, or a fixed value if specified
  !! in the global variables module.
  !! @note    The temperature should be normalised when calling this routine,
  !!          the parallel conduction is normalised on exit.
  subroutine get_kappa_para(T0_eq, tc_para)
    use mod_global_variables, only: use_fixed_tc_para, fixed_tc_para_value

    !> equilibrium temperature
    real(dp), intent(in)  :: T0_eq(gauss_gridpts)
    !> parallel thermal conduction values, normalised on exit
    real(dp), intent(out) :: tc_para(gauss_gridpts)
    real(dp)              :: T0_dimfull(gauss_gridpts)

    if (use_fixed_tc_para) then
      tc_para = fixed_tc_para_value
      return
    end if

    T0_dimfull = T0_eq * unit_temperature
    tc_para = pf_kappa_para * T0_dimfull**2.5d0 / coulomb_log
    tc_para = tc_para / unit_conduction
  end subroutine get_kappa_para


  !> Calculates the perpendicular thermal conduction.
  !! Returns either the full perpendicular thermal conduction based on
  !! the equilibrium parameters, or a fixed value if specified
  !! in the global variables module.
  !! @note    All variables should be normalised when calling this routine,
  !!          the perpendicular conduction is normalised on exit.
  subroutine get_kappa_perp(T0_eq, rho0_eq, B0_eq, tc_perp)
    use mod_global_variables, only: use_fixed_tc_perp, fixed_tc_perp_value

    !> equilibrium temperature
    real(dp), intent(in)  :: T0_eq(gauss_gridpts)
    !> equilibrium density
    real(dp), intent(in)  :: rho0_eq(gauss_gridpts)
    !> equilibrium magnetic field
    real(dp), intent(in)  :: B0_eq(gauss_gridpts)
    !> perpendicular thermal conduction values, normalised on exit
    real(dp), intent(out) :: tc_perp(gauss_gridpts)
    real(dp)              :: T0_dimfull(gauss_gridpts), B0_dimfull(gauss_gridpts)
    real(dp)              :: nH_dimfull(gauss_gridpts)

    if (use_fixed_tc_perp) then
      tc_perp = fixed_tc_perp_value
      return
    end if

    nH_dimfull = rho0_eq * unit_numberdensity
    T0_dimfull = T0_eq * unit_temperature
    B0_dimfull = B0_eq * unit_magneticfield

    tc_perp = pf_kappa_para * pf_kappa_perp * coulomb_log * nH_dimfull**2 &
              / (B0_dimfull**2 * sqrt(T0_dimfull))
    tc_perp = tc_perp / unit_conduction
  end subroutine get_kappa_perp


  !> Calculates the temperature derivative of the parallel thermal conduction component.
  !! @note    All variables should be normalised when calling this routine,
  !!          values are normalised on exit.
  subroutine get_kappa_para_derivatives(T0_eq, d_tcpara_dT)
    !> equilibrium temperature
    real(dp), intent(in)  :: T0_eq(:)
    !> temperature derivative of parallel thermal conduction, normalised on exit
    real(dp), intent(out) :: d_tcpara_dT(size(T0_eq))

    real(dp)  :: T0_dimfull(size(T0_eq))

    T0_dimfull = T0_eq * unit_temperature
    d_tcpara_dT = pf_kappa_para * 2.5d0 * T0_dimfull**1.5d0 / coulomb_log
    d_tcpara_dT = d_tcpara_dT / unit_conduction / unit_temperature
  end subroutine get_kappa_para_derivatives


  !> Calculates the thermal conduction derivatives.
  !! Returns the derivative of the perpendicular thermal conduction with respect to
  !! density, magnetic field squared and temperature.
  !! @note    All variables should be normalised when calling this routine.
  !!          The thermal conduction derivatives are normalised on exit.
  subroutine get_kappa_perp_derivatives( &
    T0_eq, rho0_eq, B0_eq, d_tc_drho, d_tc_dB2, d_tc_dT &
  )
    use mod_units, only: unit_dtc_drho, unit_dtc_dB2, unit_dtc_dT

    !> equilibrium temperature
    real(dp), intent(in)  :: T0_eq(gauss_gridpts)
    !> equilibrium density
    real(dp), intent(in)  :: rho0_eq(gauss_gridpts)
    !> equilibrium magnetic field
    real(dp), intent(in)  :: B0_eq(gauss_gridpts)
    !> derivative of perpendicular conduction with respect to density
    real(dp), intent(out) :: d_tc_drho(gauss_gridpts)
    !> derivative of perpendicular conduction with respect to magnetic field squared
    real(dp), intent(out) :: d_tc_dB2(gauss_gridpts)
    !> derivative of perpendicular conduction with respect to temperature
    real(dp), intent(out) :: d_tc_dT(gauss_gridpts)
    real(dp)    :: T0_dimfull(gauss_gridpts), rho0_dimfull(gauss_gridpts), B0_dimfull(gauss_gridpts)
    real(dp)    :: nH_dimfull(gauss_gridpts)

    ! re-dimensionalise variables, in normalised units nH = rho
    nH_dimfull = rho0_eq * unit_numberdensity
    T0_dimfull = T0_eq * unit_temperature
    rho0_dimfull = rho0_eq * unit_density
    B0_dimfull = B0_eq * unit_magneticfield

    ! density derivative
    d_tc_drho = 2.0d0 * pf_kappa_para * pf_kappa_perp * coulomb_log * nH_dimfull &
                / (B0_dimfull**2 * sqrt(T0_dimfull))
    ! magnetic field derivative
    d_tc_dB2 = -pf_kappa_para * pf_kappa_perp * coulomb_log * nH_dimfull**2 &
               / (B0_dimfull**4 * sqrt(T0_dimfull))
    ! temperature derivative
    d_tc_dT = -0.5d0 * pf_kappa_para * pf_kappa_perp * coulomb_log * nH_dimfull**2 &
              / (B0_dimfull**2 * T0_dimfull**(3.0d0/2.0d0))
    d_tc_drho = d_tc_drho / unit_dtc_drho
    d_tc_dB2 = d_tc_dB2 / unit_dtc_dB2
    d_tc_dT = d_tc_dT / unit_dtc_dT
  end subroutine get_kappa_perp_derivatives

end module mod_thermal_conduction
