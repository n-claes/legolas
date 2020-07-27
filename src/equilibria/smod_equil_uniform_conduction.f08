!
! SUBMODULE: smod_equil_uniform_conduction
!
! DESCRIPTION:
! Submodule defining an equilibrium for a uniform case, with finite thermal conduction
! in cylindrical geometry.
! Obtained from Kerner, J. Comput. Phys. 85, 1-85 (1989), Fig. 14.18
submodule (mod_equilibrium) smod_equil_uniform_conduction
  implicit none

contains

  module subroutine uniform_thermal_cond_eq()
    use mod_global_variables, only: use_fixed_tc_para, use_fixed_tc_perp, fixed_tc_para_value, fixed_tc_perp_value
    use mod_equilibrium_params, only: beta

    call allow_geometry_override(default_geometry='cylindrical', default_x_start=0.0d0, default_x_end=1.0d0)
    call initialise_grid()

    thermal_conduction = .true.
    use_fixed_tc_para = .true.
    fixed_tc_para_value = 0.001d0
    use_fixed_tc_perp = .true.
    fixed_tc_perp_value = 0.0d0

    if (use_defaults) then
      k2 = 1.0d0
      k3 = dpi

      beta = 0.25d0
    end if

    !! Parameters
    rho_field % rho0 = 1.0d0
    B_field % B02    = 0.0d0
    B_field % B03    = 1.0d0
    B_field % B0     = sqrt((B_field % B02)**2 + (B_field % B03)**2)
    ! n=1, kB=1, mu0=1 for T0
    T_field % T0     = beta * (B_field % B0)**2 / 2.0d0

  end subroutine uniform_thermal_cond_eq

end submodule smod_equil_uniform_conduction
