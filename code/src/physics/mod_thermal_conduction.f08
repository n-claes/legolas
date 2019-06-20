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
  use mod_physical_constants
  use mod_global_variables
  implicit none

  public

contains

  !> Calculates the parallel conduction coefficient as a function of the
  !! equilibrium temperature.
  !! @param[in]  T0_eq    Array containing the equilibrium temperatures, in K
  !! @param[out] tc_para  Array containing the parallel conduction coefficients
  subroutine get_kappa_para(T0_eq, tc_para)
    real(dp), intent(in)  :: T0_eq(4*gridpts)
    real(dp), intent(out) :: tc_para(4*gridpts)

    real(dp)              :: prefactor_para

    if (cgs_units) then
      prefactor_para = 1.8d-5
    else
      prefactor_para = 1.8d-10
    end if

    tc_para = prefactor_para * T0_eq**2.5 / coulomb_log

  end subroutine get_kappa_para

  !> Calculates the perpendicular conduction coefficient as a function of the
  !! equilibrium temperature, density and magnetic field.
  !! @param[in]  T0_eq    Array containing the equilibrium temperatures, in K
  !! @param[in]  rho0_eq  Array containing the equilibrium density
  !! @param[in]  B0_eq    Array containing the equilibrium magnetic field
  !! @param[out] tc_perp  Array containing the perpendicular conduction coefficients
  subroutine get_kappa_perp(T0_eq, rho0_eq, B0_eq, tc_perp)
    real(dp), intent(in)  :: T0_eq(4*gridpts), rho0_eq(4*gridpts), &
                             B0_eq(4*gridpts)
    real(dp), intent(out) :: tc_perp(4*gridpts)

    real(dp)              :: tc_para(4*gridpts), nH(4*gridpts)
    real(dp)              :: prefactor_perp, mp

    if (cgs_units) then
      mp = mp_cgs
    else
      mp = mp_si
    end if

    !! cgs is already accounted for when calling get_kappa_para
    prefactor_perp = 8.2d-33

    call get_kappa_para(T0_eq, tc_para)
    call get_nH(rho0_eq, mp, nH)
    tc_perp = prefactor_perp * coulomb_log**2 * nH**2 * tc_para &
              / (B0_eq**2 * T0_eq**3)

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
    real(dp), intent(in)  :: T0_eq(4*gridpts), rho0_eq(4*gridpts), &
                             B0_eq(4*gridpts)
    real(dp), intent(out) :: d_tc_drho(4*gridpts)

    real(dp)              :: nH(4*gridpts)
    real(dp)              :: prefactor_perp, mp

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
    d_tc_drho = 2*prefactor_perp * coulomb_log * nH / (B0_eq**2 * T0_eq**0.5)

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
    real(dp), intent(in)  :: T0_eq(4*gridpts), rho0_eq(4*gridpts), &
                             B0_eq(4*gridpts)
    real(dp), intent(out) :: d_tc_dT(4*gridpts)

    real(dp)              :: nH(4*gridpts)
    real(dp)              :: prefactor_perp, mp

    if (cgs_units) then
      prefactor_perp = 1.476d-37
      mp = mp_cgs
    else
      prefactor_perp = 1.476d-42
      mp = mp_si
    end if

    call get_nH(rho0_eq, mp, nH)
    d_tc_dT = -0.5d0*prefactor_perp * coulomb_log * nH**2 &
              / (B0_eq**2 * T0_eq**1.5)

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
    real(dp), intent(in)  :: T0_eq(4*gridpts), rho0_eq(4*gridpts), &
                             B0_eq(4*gridpts)
    real(dp), intent(out) :: d_tc_dB2(4*gridpts)

    real(dp)              :: nH(4*gridpts)
    real(dp)              :: prefactor_perp, mp

    if (cgs_units) then
      prefactor_perp = 1.476d-37
      mp = mp_cgs
    else
      prefactor_perp = 1.476d-42
      mp = mp_si
    end if

    call get_nH(rho0_eq, mp, nH)
    d_tc_dB2 = -1*prefactor_perp * coulomb_log * nH**2 &
              / (B0_eq**4 * T0_eq**0.5)

  end subroutine get_dkappa_perp_dB2

  !> Calculates the hydrogen number density based on the equilibrium density.
  !! @param[in]  rho0_eq   Array containing the equilibrium density
  !! @param[in]  mp        The proton mass, in cgs or SI units
  !! @param[out] nH        Array containing the equilibrium number density
  subroutine get_nH(rho0_eq, mp, nH)
    real(dp), intent(in)  :: rho0_eq(4*gridpts), mp
    real(dp), intent(out) :: nH(4*gridpts)

    nH = rho0_eq / ((1.0d0 + 4.0d0 * He_abundance) * mp)
  end subroutine get_nH


end module mod_thermal_conduction
