!
! SUBMODULE: smod_equil_magnetothermal_instabilities
!
! DESCRIPTION:
!> Submodule defining magnetothermal instabilities in a cylindrical geometry as they appear in
!! Van Der Linden et al., Solar Physics 140, 317-342, 1992. Uses the radiative cooling curve defined by
!! Rosner et al. (1978) and includes (parallel) thermal conduction.
submodule(mod_equilibrium) smod_equil_magnetothermal_instabilities
  implicit none

contains

  module subroutine magnetothermal_instability_eq()
    use mod_equilibrium_params, only: cte_T0
    use mod_global_variables, only: cooling_curve, use_fixed_tc_perp, fixed_tc_perp_value

    real(dp)  :: r, p_r(gauss_gridpts)
    integer   :: i

    geometry = 'cylindrical'
    x_start = 0.0d0
    x_end = 1.0d0
    call initialise_grid()

    radiative_cooling = .true.
    cooling_curve = 'rosner'
    thermal_conduction = .true.
    use_fixed_tc_perp = .true.
    fixed_tc_perp_value = 0.0d0

    if (use_defaults) then
      cte_T0 = 1.0d0

      k2 = 0.0d0
      k3 = 1.0d0
    end if

    cgs_units = .true.
    call define_temp_mag_len(new_unit_temperature = 2.6d6, &
                             new_unit_magneticfield = 10.0d0, &
                             new_unit_length = 2.44d9)

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      B_field % B02(i) = r / (1.0d0 + r**2)
      B_field % B03(i) = 0.0d0
      B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      p_r(i) = 1.0d0 / ( 2.0d0 * (1.0d0 + r**2)**2 )
      T_field % T0(i) = cte_T0
      rho_field % rho0(i) = p_r(i) / cte_T0

      rho_field % d_rho0_dr(i) = -2.0d0 * r / (cte_T0 * (r**2 + 1.0d0)**3)
      B_field % d_B02_dr(i) = (1.0d0 - r**2) / (r**4 + 2.0d0 * r**2 + 1.0d0)
      ! Temperature derivative is zero due to isothermal
    end do

  end subroutine magnetothermal_instability_eq

end submodule smod_equil_magnetothermal_instabilities
