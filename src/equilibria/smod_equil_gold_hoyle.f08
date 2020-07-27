!
! SUBMODULE: smod_equil_gold_hoyle
!
! DESCRIPTION:
! Submodule defining a Gold-Hoyle equilibrium in cylindrical geometry. This equilibrium configuration
! models a filament with a uniform twist, such that all fieldlines perform an equal amount of turns around
! the cylinder axis. Obtained from Van Der Linden, Goossens and Hood, SoPh, 140, 317V (1992).
! Three profiles:
! - cold: B = 10.0 G, rho = 1.6726e-12 kg/m3, R = 1e8 m, T = 5.7e4 K
! - cool: B = 22.5 G, rho = 1.6726e-12 kg/m3, R = 1e8 m, T = 2.9e5 K
! - hot : B = 67.0 G, rho = 1.6726e-12 kg/m3, R = 1e9 m, T = 2.6e6 K
submodule (mod_equilibrium) smod_equil_gold_hoyle
  implicit none

contains

  module subroutine gold_hoyle_eq()
    use mod_equilibrium_params, only: cte_T0, cte_rho0, alpha
    use mod_global_variables, only: cooling_curve, use_fixed_tc_perp, fixed_tc_perp_value

    real(dp)  :: r
    integer   :: i

    call allow_geometry_override(default_geometry='cylindrical', default_x_start=0.0d0, default_x_end=1.0d0)
    call initialise_grid()

    if (use_defaults) then
      ! physics
      radiative_cooling = .true.
      cooling_curve = "rosner"
      thermal_conduction = .true.
      use_fixed_tc_perp = .true.
      fixed_tc_perp_value = 0.0d0

      ! parameters
      k2 = 1.0d0
      k3 = 1.0d0
      cte_rho0 = 1.0d0
      cte_T0 = 0.001d0
      alpha = 20.0d0

      ! units
      cgs_units = .true.
      ! cool plasma
      call set_normalisations(new_unit_density=1.6727d-15, new_unit_magneticfield=22.5d0, new_unit_length=1.0d10)
      ! hot plasma
      ! call set_normalisations(new_unit_density=1.6727d-15, new_unit_magneticfield=67.0d0, new_unit_length=1.0d11)
      ! cold plasma
      ! call set_normalisations(new_unit_density=1.6727d-15, new_unit_magneticfield=10.0d0, new_unit_length=1.0d10)
    end if

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      rho_field % rho0(i) = cte_rho0
      T_field % T0(i) = cte_T0
      B_field % B02(i) = alpha * r / (1.0d0 + alpha**2 * r**2)
      B_field % B03(i) = 1.0d0 / (1.0d0 + alpha**2 * r**2)
      B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)

      B_field % d_B02_dr(i) = -(alpha**3 * r**2 - alpha) / (alpha**4 * r**4 + 2.0d0 * alpha**2 * r**2 + 1.0d0)
      B_field % d_B03_dr(i) = -2.0d0 * alpha**2 * r / (alpha**2 * r**2 + 1.0d0)**2
    end do

  end subroutine gold_hoyle_eq

end submodule smod_equil_gold_hoyle
