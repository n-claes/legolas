!
! MODULE: mod_thermal_conduction
!
!> @author
!> Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Module to calculate the thermal conduction contributions.
!
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

  subroutine set_conduction_values(rho_field, T_field, B_field, kappa_field)
    use mod_types, only: density_type, temperature_type, bfield_type, conduction_type
    use mod_global_variables, only: use_fixed_tc_perp

    type(density_type), intent(in)        :: rho_field
    type(temperature_type), intent(in)    :: T_field
    type(bfield_type), intent(in)         :: B_field
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
    call get_kappa_perp(T_field % T0, rho_field % rho0, B_field % B0, kappa_field % kappa_perp)

    if (.not. use_fixed_tc_perp) then
      call get_kappa_perp_derivatives(T_field % T0, rho_field % rho0, B_field % B0, &
                                      kappa_field % d_kappa_perp_drho, &
                                      kappa_field % d_kappa_perp_dB2, &
                                      kappa_field % d_kappa_perp_dT)
    end if
  end subroutine set_conduction_values


  !> Calculates the parallel conduction coefficient as a function of the
  !! equilibrium temperature.
  !! @param[in]  T0_eq    normalised equilibrium temperature
  !! @param[out] tc_para  normalised parallel thermal conduction
  subroutine get_kappa_para(T0_eq, tc_para)
    use mod_global_variables, only: use_fixed_tc_para, fixed_tc_para_value

    real(dp), intent(in)  :: T0_eq(gauss_gridpts)
    real(dp), intent(out) :: tc_para(gauss_gridpts)
    real(dp)              :: T0_dimfull(gauss_gridpts)

    if (use_fixed_tc_para) then
      tc_para = fixed_tc_para_value
      return
    end if

    T0_dimfull = T0_eq * unit_temperature
    tc_para = pf_kappa_para * T0_dimfull**2.5 / coulomb_log
    tc_para = tc_para / unit_conduction
  end subroutine get_kappa_para

  !> Calculates the perpendicular conduction coefficient as a function of the
  !! equilibrium temperature, density and magnetic field.
  !! @param[in]  T0_eq    normalised equilibrium temperature
  !! @param[in]  rho0_eq  normalised equilibrium density
  !! @param[in]  B0_eq    normalised equilibrium magnetic field
  !! @param[out] tc_perp  normalised perpendicular thermal conduction
  subroutine get_kappa_perp(T0_eq, rho0_eq, B0_eq, tc_perp)
    use mod_global_variables, only: use_fixed_tc_perp, fixed_tc_perp_value

    real(dp), intent(in)  :: T0_eq(gauss_gridpts), rho0_eq(gauss_gridpts), B0_eq(gauss_gridpts)
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


  subroutine get_kappa_perp_derivatives(T0_eq, rho0_eq, B0_eq, d_tc_drho, d_tc_dB2, d_tc_dT)
    use mod_units, only: unit_dtc_drho, unit_dtc_dB2, unit_dtc_dT

    real(dp), intent(in)  :: T0_eq(gauss_gridpts), rho0_eq(gauss_gridpts), B0_eq(gauss_gridpts)
    real(dp), intent(out) :: d_tc_drho(gauss_gridpts), d_tc_dB2(gauss_gridpts), d_tc_dT(gauss_gridpts)
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
              / (B0_dimfull**2 * T0_dimfull**(-3.0d0/2.0d0))
    d_tc_drho = d_tc_drho / unit_dtc_drho
    d_tc_dB2 = d_tc_dB2 / unit_dtc_dB2
    d_tc_dT = d_tc_dT / unit_dtc_dT
  end subroutine get_kappa_perp_derivatives

end module mod_thermal_conduction
