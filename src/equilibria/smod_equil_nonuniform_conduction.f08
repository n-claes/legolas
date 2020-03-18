!
! SUBMODULE: smod_equil_nonuniform_conduction
!
! DESCRIPTION:
!> Submodule defining an equilibrium for a non-uniform case, with finite thermal_conduction
!! conduction in cylindrical geometry.
!! Obtained from Kerner, J. Comput. Phys. 85, 1-85 (1989), Fig. 14.19
submodule (mod_equilibrium) smod_equil_nonuniform_conduction
  implicit none

contains

  module subroutine nonuniform_thermal_cond_eq()
    use mod_global_variables, only: use_fixed_tc_para, use_fixed_tc_perp, fixed_tc_para_value, fixed_tc_perp_value
    use mod_equilibrium_params, only: beta

    real(dp)      :: x, B0
    integer       :: i

    geometry = 'Cartesian'
    x_start = 0.0d0
    x_end   = 1.0d0
    call initialise_grid()

    thermal_conduction = .true.
    use_fixed_tc_para = .true.
    fixed_tc_para_value = 0.001d0
    use_fixed_tc_perp = .true.
    fixed_tc_perp_value = 0.0d0

    if (use_defaults) then
      k2 = 0.0d0
      k3 = 1.0d0
      beta = 0.25d0
    end if

    B0 = 1.0d0

    rho_field % rho0 = 1.0d0
    T_field % T0 = 1.0d0

    do i = 1, gauss_gridpts
      x = grid_gauss(i)

      !! Equilibrium
      B_field % B02(i) = B0*sin(x)
      B_field % B03(i) = B0*cos(x)
      B_field % B0(i)  = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)

      !! Derivatives
      B_field % d_B02_dr(i) = B0*cos(x)
      B_field % d_B03_dr(i) = -B0*sin(x)
    end do

  end subroutine nonuniform_thermal_cond_eq

end submodule smod_equil_nonuniform_conduction
