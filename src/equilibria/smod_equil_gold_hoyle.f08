! =============================================================================
!> This submodule defines a Gold-Hoyle equilibrium in cylindrical geometry.
!! This equilibrium configuration models a filament with a uniform twist
!! such that all fieldlines perform an equal amount of turns around the cylinder axis.
!! The geometry can be overridden in the parfile.
!!
!! This equilibrium is taken from
!! _Van der Linden, R. A. M., Goossens, M., & Hood, A. W. (1992).
!!  The relevance of the ballooning approximation for magnetic, thermal,
!!  and coalesced magnetothermal instabilities. Solar physics, 140(2), 317-342_.
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 1
!! - <tt>k3</tt> = 1
!! - <tt>cte_T0</tt> = 1 : used to set the temperature.
!! - <tt>cte_rho0</tt> = 1 : used to set the density
!! - <tt>alpha</tt> = 20 : used in the magnetic field components
!! - cooling_curve = 'rosner'
!! - parallel thermal conduction, no perpendicular conduction
!!
!! and can all be changed in the parfile. @endnote
!! @note To reproduce the three profiles of the original paper you can supply
!!       one of the following normalisations:
!!
!! - cold plasma: <tt>unit_magneticfield</tt> = 10 Gauss,
!!                <tt>unit_density</tt> = 1.6726e-15 g/cm3,
!!                <tt>unit_length</tt> = 1.0e10 cm,
!!                corresponding to a temperature of 5.7e4 K.
!! - cool plasma: <tt>unit_magneticfield</tt> = 22.5 Gauss,
!!                <tt>unit_density</tt> = 1.6726e-15 g/cm3,
!!                <tt>unit_length</tt> = 1.0e10 cm,
!!                corresponding to a temperature of 2.9e5 K (default).
!! - hot plasma:  <tt>unit_magneticfield</tt> = 67 Gauss,
!!                <tt>unit_density</tt> = 1.6726e-15 g/cm3,
!!                <tt>unit_length</tt> = 1.0e11 cm,
!!                corresponding to a temperature of 2.6e6 K. @endnote
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
      radiative_cooling = .true.
      cooling_curve = "rosner"
      thermal_conduction = .true.
      use_fixed_tc_perp = .true.
      fixed_tc_perp_value = 0.0d0

      k2 = 1.0d0
      k3 = 1.0d0
      cte_rho0 = 1.0d0
      cte_T0 = 0.001d0
      alpha = 20.0d0

      cgs_units = .true.
      ! cool plasma: B = 22.5 G, rho = 1.6726e-12 kg/m3, R = 1e8 m, T = 2.9e5 K
      call set_normalisations(new_unit_density=1.6727d-15, new_unit_magneticfield=22.5d0, new_unit_length=1.0d10)
      ! hot plasma: B = 67.0 G, rho = 1.6726e-12 kg/m3, R = 1e9 m, T = 2.6e6 K
      ! call set_normalisations(new_unit_density=1.6727d-15, new_unit_magneticfield=67.0d0, new_unit_length=1.0d11)
      ! cold plasma: B = 10.0 G, rho = 1.6726e-12 kg/m3, R = 1e8 m, T = 5.7e4 K
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
