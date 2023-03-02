! =============================================================================
!> This submodule defines a resistive equilibrium with tearing modes created
!! by a Harris sheet.
!!
!! This equilibrium is taken from Shi et al. (2020), Oblique tearing mode
!! instability: guide field and Hall effect. _The Astrophysical Journal_, 902:142.
!! [DOI](https://doi.org/10.3847/1538-4357/abb6fa).
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 0.12
!! - <tt>k3</tt> = 0
!! - <tt>cte_rho0</tt> = 1
!! - <tt>cte_T0</tt> = 1
!! - <tt>cte_B02</tt> = 1
!! - <tt>cte_B03</tt> = 0 : guide field parameter.
!! - <tt>alpha</tt> = 1 : used to set the width of the current sheet.
!! - <tt>eq_bool</tt> = False : if True, an alternative force-free Harris sheet is used.
!!
!! and can all be changed in the parfile. @endnote
submodule (mod_equilibrium) smod_equil_harris_sheet
  implicit none

contains

  module procedure harris_sheet_eq
    use mod_equilibrium_params, only: cte_rho0, cte_B02, cte_B03, cte_T0, alpha, eq_bool

    real(dp)      :: x
    integer       :: i

    if (settings%equilibrium%use_defaults) then
      call settings%grid%set_geometry("Cartesian")
      call settings%grid%set_grid_boundaries(-15.0_dp, 15.0_dp)
      call settings%physics%enable_resistivity(fixed_resistivity_value=0.001_dp)

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
    call initialise_grid(settings)

    do i = 1, settings%grid%get_gauss_gridpts()
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
      do i = 1, settings%grid%get_gauss_gridpts()
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
      do i = 1, settings%grid%get_gauss_gridpts()
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
  end procedure harris_sheet_eq

end submodule smod_equil_harris_sheet
