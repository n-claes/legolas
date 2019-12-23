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
  use mod_physical_constants, only: dpi, coulomb_log, mp_cgs, mp_si, unit_temperature
  implicit none

  private
  
  public :: get_kappa_para
  public :: get_kappa_perp
  public :: get_dkappa_perp_drho
  public :: get_dkappa_perp_dT
  public :: get_dkappa_perp_dB2

contains

  !> Calculates the parallel conduction coefficient as a function of the
  !! equilibrium temperature.
  !! @param[in]  T0_eq    Array containing the equilibrium temperatures, in K
  !! @param[out] tc_para  Array containing the parallel conduction coefficients
  subroutine get_kappa_para(T0_eq, tc_para)
    use mod_physical_constants, only: unit_conduction
    
    real(dp), intent(in)  :: T0_eq(gauss_gridpts)
    real(dp), intent(out) :: tc_para(gauss_gridpts)

    real(dp)              :: prefactor_para
    real(dp)              :: T0_eq_denorm(gauss_gridpts)

    !! Denormalise variables for calculation
    T0_eq_denorm = T0_eq * unit_temperature

    if (cgs_units) then
      prefactor_para = 1.8d-5
    else
      prefactor_para = 1.8d-10
    end if

    tc_para = prefactor_para * T0_eq_denorm**2.5 / coulomb_log
    !! Renormalise
    tc_para = tc_para / unit_conduction

  end subroutine get_kappa_para

  !> Calculates the perpendicular conduction coefficient as a function of the
  !! equilibrium temperature, density and magnetic field.
  !! @param[in]  T0_eq    Array containing the equilibrium temperatures, in K
  !! @param[in]  rho0_eq  Array containing the equilibrium density
  !! @param[in]  B0_eq    Array containing the equilibrium magnetic field
  !! @param[out] tc_perp  Array containing the perpendicular conduction coefficients
  subroutine get_kappa_perp(T0_eq, rho0_eq, B0_eq, tc_perp)
    use mod_physical_constants, only: unit_magneticfield, unit_conduction, unit_numberdensity
    
    real(dp), intent(in)  :: T0_eq(gauss_gridpts), rho0_eq(gauss_gridpts), &
                             B0_eq(gauss_gridpts)
    real(dp), intent(out) :: tc_perp(gauss_gridpts)

    real(dp)              :: T0_eq_denorm(gauss_gridpts), &
                             B0_eq_denorm(gauss_gridpts)
    real(dp)              :: tc_para(gauss_gridpts), nH(gauss_gridpts)
    real(dp)              :: prefactor_perp, mp

    !! Denormalise variables for calculation
    T0_eq_denorm   = T0_eq * unit_temperature
    B0_eq_denorm   = B0_eq * unit_magneticfield

    if (cgs_units) then
      mp = mp_cgs
    else
      mp = mp_si
    end if

    !! cgs is already accounted for when calling get_kappa_para
    prefactor_perp = 8.2d-33

    call get_kappa_para(T0_eq, tc_para)
    tc_para = tc_para * unit_conduction

    call get_nH(rho0_eq, mp, nH)
    nH = nH * unit_numberdensity

    tc_perp = prefactor_perp * coulomb_log**2 * nH**2 * tc_para &
              / (B0_eq_denorm**2 * T0_eq_denorm**3)
    !! Renormalise
    tc_perp = tc_perp / unit_conduction

  end subroutine get_kappa_perp

  !> Calculates the derivative of the perpendicular conduction coefficient
  !! with respect to density.
  !! @param[in]  T0_eq    Array containing the equilibrium temperatures, in K
  !! @param[in]  rho0_eq  Array containing the equilibrium density
  !! @param[in]  B0_eq    Array containing the equilibrium magnetic field
  !! @param[out] d_tc_drho  Array containing the derivative of the
  !!                        perpendicular conduction coefficient
  !!                        with respect to density
  subroutine get_dkappa_perp_drho(T0_eq, rho0_eq, B0_eq, d_tc_drho)
    use mod_physical_constants, only: unit_magneticfield, unit_dtc_drho, unit_numberdensity
    
    real(dp), intent(in)  :: T0_eq(gauss_gridpts), rho0_eq(gauss_gridpts), &
                             B0_eq(gauss_gridpts)
    real(dp), intent(out) :: d_tc_drho(gauss_gridpts)

    real(dp)              :: T0_eq_denorm(gauss_gridpts), &
                             B0_eq_denorm(gauss_gridpts)
    real(dp)              :: nH(gauss_gridpts)
    real(dp)              :: prefactor_perp, mp

    !! Denormalise variables for calculation
    T0_eq_denorm = T0_eq * unit_temperature
    B0_eq_denorm = B0_eq * unit_magneticfield

    if (cgs_units) then
      !! 1.8d-10 * 8.2d-33 * 1.0d5
      prefactor_perp = 1.476d-37
      mp = mp_cgs
    else
      !! 1.8d-10 * 8.2d-33
      prefactor_perp = 1.476d-42
      mp = mp_si
    end if

    call get_nH(rho0_eq, mp, nH)
    nH = nH * unit_numberdensity

    d_tc_drho = 2*prefactor_perp * coulomb_log * nH / &
                      (B0_eq_denorm**2 * T0_eq_denorm**0.5)
    !! Renormalise
    d_tc_drho = d_tc_drho / unit_dtc_drho

  end subroutine get_dkappa_perp_drho

  !> Calculates the derivative of the perpendicular conduction coefficient
  !! with respect to temperature.
  !! @param[in]  T0_eq    Array containing the equilibrium temperatures, in K
  !! @param[in]  rho0_eq  Array containing the equilibrium density
  !! @param[in]  B0_eq    Array containing the equilibrium magnetic field
  !! @param[out] d_tc_dT  Array containing the derivative of the
  !!                      perpendicular conduction coefficient
  !!                      with respect to temperature
  subroutine get_dkappa_perp_dT(T0_eq, rho0_eq, B0_eq, d_tc_dT)
    use mod_physical_constants, only: unit_magneticfield, unit_dtc_dT, unit_numberdensity
    
    real(dp), intent(in)  :: T0_eq(gauss_gridpts), rho0_eq(gauss_gridpts), &
                             B0_eq(gauss_gridpts)
    real(dp), intent(out) :: d_tc_dT(gauss_gridpts)

    real(dp)              :: T0_eq_denorm(gauss_gridpts), &
                             B0_eq_denorm(gauss_gridpts)
    real(dp)              :: nH(gauss_gridpts)
    real(dp)              :: prefactor_perp, mp

    !! Denormalise variables for calculation
    T0_eq_denorm = T0_eq * unit_temperature
    B0_eq_denorm = B0_eq * unit_magneticfield

    if (cgs_units) then
      prefactor_perp = 1.476d-37
      mp = mp_cgs
    else
      prefactor_perp = 1.476d-42
      mp = mp_si
    end if

    call get_nH(rho0_eq, mp, nH)
    nH = nH * unit_numberdensity

    d_tc_dT = -0.5d0*prefactor_perp * coulomb_log * nH**2 &
              / (B0_eq_denorm**2 * T0_eq_denorm**1.5)
    !! Renormalise
    d_tc_dT = d_tc_dT / unit_dtc_dT

  end subroutine get_dkappa_perp_dT

  !> Calculates the derivative of the perpendicular conduction coefficient
  !! with respect to B^2.
  !! @param[in]  T0_eq    Array containing the equilibrium temperatures, in K
  !! @param[in]  rho0_eq  Array containing the equilibrium density
  !! @param[in]  B0_eq    Array containing the equilibrium magnetic field
  !! @param[out] d_tc_dB2  Array containing the derivative of the
  !!                       perpendicular conduction coefficient
  !!                       with respect to B^2
  subroutine get_dkappa_perp_dB2(T0_eq, rho0_eq, B0_eq, d_tc_dB2)
    use mod_physical_constants, only: unit_magneticfield, unit_dtc_dB2, unit_numberdensity
    
    real(dp), intent(in)  :: T0_eq(gauss_gridpts), rho0_eq(gauss_gridpts), &
                             B0_eq(gauss_gridpts)
    real(dp), intent(out) :: d_tc_dB2(gauss_gridpts)

    real(dp)              :: T0_eq_denorm(gauss_gridpts), &
                             B0_eq_denorm(gauss_gridpts)
    real(dp)              :: nH(gauss_gridpts)
    real(dp)              :: prefactor_perp, mp

    !! Denormalise variables for calculation
    T0_eq_denorm = T0_eq * unit_temperature
    B0_eq_denorm = B0_eq * unit_magneticfield

    if (cgs_units) then
      prefactor_perp = 1.476d-37
      mp = mp_cgs
    else
      prefactor_perp = 1.476d-42
      mp = mp_si
    end if

    call get_nH(rho0_eq, mp, nH)
    nH = nH * unit_numberdensity

    d_tc_dB2 = -1*prefactor_perp * coulomb_log * nH**2 &
              / (B0_eq_denorm**4 * T0_eq_denorm**0.5)
    !! Renormalise
    d_tc_dB2 = d_tc_dB2 / unit_dtc_dB2

  end subroutine get_dkappa_perp_dB2

  !> Calculates the hydrogen number density based on the equilibrium density.
  !! @param[in]  rho0_eq   Array containing the equilibrium density
  !! @param[in]  mp        The proton mass, in cgs or SI units
  !! @param[out] nH        Array containing the equilibrium number density
  subroutine get_nH(rho0_eq, mp, nH)
    use mod_physical_constants, only: unit_density, unit_numberdensity, He_abundance
    
    real(dp), intent(in)  :: rho0_eq(gauss_gridpts), mp
    real(dp), intent(out) :: nH(gauss_gridpts)
    real(dp)              :: rho0_eq_denorm(gauss_gridpts)

    !! Denormalise variables for calculation
    rho0_eq_denorm = rho0_eq * unit_density

    nH = rho0_eq_denorm / ((1.0d0 + 4.0d0 * He_abundance) * mp)
    !! Renormalise
    nH = nH / unit_numberdensity
  end subroutine get_nH


end module mod_thermal_conduction
