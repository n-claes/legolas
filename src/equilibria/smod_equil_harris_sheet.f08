! =============================================================================
!> This submodule defines a resistive equilibrium with tearing modes created
!! by a Harris sheet.
!!
!! This equilibrium is taken from Shi et al. (2020), Oblique tearing mode
!! instability: guide field and Hall effect
submodule (mod_equilibrium) smod_equil_harris_sheet
  implicit none

contains

  module subroutine harris_sheet_eq()
    use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value
    use mod_equilibrium_params, only: cte_rho0, cte_B02, cte_B03, cte_T0, alpha, eq_bool

    real(dp)      :: x
    integer       :: i

    call allow_geometry_override(default_geometry='Cartesian', default_x_start=-15.0d0, default_x_end=15.0d0)
    call initialise_grid()

    if (use_defaults) then
      resistivity = .true.
      use_fixed_resistivity = .true.
      fixed_eta_value = 0.0001d0

      k2 = 0.12d0
      k3 = 0.0d0

      alpha = 1.0d0

      cte_rho0 = 1.0d0
      cte_B02 = 1.0d0
      cte_B03 = 0.0d0
      cte_T0 = 1.0d0

      !> eq_bool >> if True, the alternative force-free Harris sheet is used
      eq_bool = .false.
    end if

    do i = 1, gauss_gridpts
      x = grid_gauss(i)

      rho_field % rho0(i) = cte_rho0
      B_field % B02(i)    = cte_B02 * tanh(x / alpha)

      rho_field % d_rho0_dr(i) = 0.0d0
      B_field % d_B02_dr(i)    = cte_B02 / (alpha * cosh(x / alpha)**2)

      eta_field % dd_B02_dr(i) = - 2.0d0 * cte_B02 * sinh(x / alpha) &
                                  / (alpha**2 * cosh(x / alpha)**3)
    end do

    if (eq_bool) then
      !> force-free option: B0 = constant, so T0 can be any arbitrary constant
      do i = 1, gauss_gridpts
        x = grid_gauss(i)

        B_field % B03(i)    = sqrt(cte_B03**2 + cte_B02**2 / cosh(x / alpha)**2)
        B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
        T_field % T0(i)     = cte_T0

        B_field % d_B03_dr(i)    = -cte_B02**2 * sinh(x / alpha) &
                                    / (alpha * cosh(x / alpha)**3 * (B_field % B03(i)))
        T_field % d_T0_dr(i)     = 0.0d0
        eta_field % dd_B03_dr(i) = cte_B02**2 * (-1.0d0 + 2.0d0 * sinh(x / alpha)**2 &
                                    - cte_B02**2 * tanh(x / alpha)**2 &
                                    / (B_field % B03(i))**2 &
                                    ) / (alpha**2 * cosh(x / alpha)**4 * (B_field % B03(i)))
      end do
    else
      do i = 1, gauss_gridpts
        x = grid_gauss(i)

        B_field % B03(i)    = cte_B03
        B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
        T_field % T0(i)     = (cte_B03**2 + cte_B02**2 - (B_field % B0(i))**2) &
                                / (2.0d0 * cte_rho0)

        B_field % d_B03_dr(i)    = 0.0d0
        T_field % d_T0_dr(i)     = -cte_B02**2 * sinh(x / alpha) &
                                    & / (alpha * cte_rho0 * cosh(x / alpha)**3)
        eta_field % dd_B03_dr(i) = 0.0d0
      end do
    end if
  end subroutine harris_sheet_eq

end submodule smod_equil_harris_sheet
