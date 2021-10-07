! =============================================================================
!> This submodule defines an equilibrium containing magnetothermal instabilities
!! in a cylindrical geometry. The geometry can be overridden in the parfile.
!!
!! This equilibrium is taken from
!! _Van der Linden, R. A. M., Goossens, M., & Hood, A. W. (1992).
!!  The relevance of the ballooning approximation for magnetic, thermal,
!!  and coalesced magnetothermal instabilities. Solar physics, 140(2), 317-342_.
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 0
!! - <tt>k3</tt> = 1
!! - <tt>cte_T0</tt> = 1 : used to set the temperature.
!! - cooling_curve = 'rosner'
!! - parallel thermal conduction, no perpendicular conduction
!!
!! and normalisations given by
!!
!! - <tt>unit_temperature</tt> = 2.6e6 K
!! - <tt>unit_magneticfield</tt> = 10 Gauss
!! - <tt>unit_length</tt> = 1e8 cm
!!
!! and can all be changed in the parfile. @endnote
!! @note The default setup handles _CASE I_ in the original paper. For _CASE II_,
!!       you can set k2 = 10. To reproduce the thermal continuum plot,
!!       set <tt>unit_length</tt> = 2.44e9 cm in _CASE I_. @endnote
submodule(mod_equilibrium) smod_equil_magnetothermal_instabilities
  implicit none

contains

  !> Sets the equilibrium
  module subroutine magnetothermal_instability_eq()
    use mod_equilibrium_params, only: cte_T0
    use mod_global_variables, only: cooling_curve, use_fixed_tc_perp, fixed_tc_perp_value

    real(dp)  :: r, p_r(gauss_gridpts)
    integer   :: i

    call allow_geometry_override( &
      default_geometry="cylindrical", default_x_start=0.0d0, default_x_end=1.0d0 &
    )
    call initialise_grid()

    if (use_defaults) then ! LCOV_EXCL_START
      radiative_cooling = .true.
      cooling_curve = "rosner"
      thermal_conduction = .true.
      use_fixed_tc_perp = .true.
      fixed_tc_perp_value = 0.0d0

      cgs_units = .true.
      call set_normalisations( &
        new_unit_temperature=2.6d6, &
        new_unit_magneticfield=10.0d0, &
        new_unit_length=1.00d8, &
        new_mean_molecular_weight=1.0d0 & ! this work assumes pure proton plasma
      )

      cte_T0 = 1.0d0
      k2 = 0.0d0
      k3 = 1.0d0
    end if ! LCOV_EXCL_STOP

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
