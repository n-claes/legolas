module mod_thermal_conduction
  use mod_global_variables
  implicit none

contains

  subroutine calculate_para_conduction(tc_para)
    use mod_physical_constants

    real(dp), intent(out) :: tc_para(4*gridpts)
    real(dp)              :: prefactor_para

    if (cgs_units) then
      prefactor_para = 8.0d-7
    else
      prefactor_para = 8.0d-12
    end if

    ! TODO
    tc_para = 1.0d0

  end subroutine calculate_para_conduction

  subroutine calculate_perp_conduction(tc_perp)
    use mod_physical_constants

    real(dp), intent(out) :: tc_perp(4*gridpts)
    real(dp)              :: prefactor_perp

    if (cgs_units) then
      prefactor_perp = 4.0d-10
    else
      prefactor_perp = 4.0d-30
    end if

    ! TODO
    tc_perp = 1.0d0
  end subroutine calculate_perp_conduction


end module mod_thermal_conduction
