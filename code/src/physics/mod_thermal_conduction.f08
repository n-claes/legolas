module mod_thermal_conduction
  use mod_global_variables
  implicit none

  real(dp), parameter       :: coulomb_log = 22

contains

  subroutine calculate_para_conduction(T0_eq, tc_para)
    use mod_physical_constants

    real(dp), intent(in)  :: T0_eq(4*gridpts)
    real(dp), intent(out) :: tc_para(4*gridpts)

    real(dp)              :: prefactor_para

    if (cgs_units) then
      prefactor_para = 1.8d-5
    else
      prefactor_para = 1.8d-10
    end if

    tc_para = prefactor_para * T0_eq**2.5 / coulomb_log

  end subroutine calculate_para_conduction

  subroutine calculate_perp_conduction(T0_eq, rho0_eq, B0_eq, tc_perp)
    use mod_physical_constants

    real(dp), intent(in)  :: T0_eq(4*gridpts), rho0_eq(4*gridpts), &
                             B0_eq(4*gridpts)
    real(dp), intent(out) :: tc_perp(4*gridpts)

    real(dp)              :: tc_para(4*gridpts), nH(4*gridpts)
    real(dp)              :: prefactor_perp, mp

    if (cgs_units) then
      prefactor_perp = 8.2d-27
      mp = mp_cgs
    else
      prefactor_perp = 8.2d-33
      mp = mp_si
    end if

    call calculate_para_conduction(T0_eq, tc_para)
    nH = rho0_eq / ((1.0d0 + 4.0d0 * He_abundance) * mp)
    tc_perp = coulomb_log**2 * nH**2 * tc_para / (B0_eq**2 * T0_eq**3)

  end subroutine calculate_perp_conduction


end module mod_thermal_conduction
