! =============================================================================
!> This submodule defines a coronal loop modelled as a cylinder with viscous
!! and resistive effects included.
!!
!! This equilibrium is taken from Erdélyi & Goossens (1995)
! SUBMODULE: smod_equil_coronal_flux_tube
submodule(mod_equilibrium) smod_equil_viscoresistive_tube
  implicit none

contains

  module subroutine viscoresistive_tube_eq()
    use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value, viscosity_value
    use mod_equilibrium_params, only: cte_B02, alpha, delta

    real(dp)  :: r
    real(dp)  :: p0(gauss_gridpts), dp0(gauss_gridpts)
    integer   :: i

    call allow_geometry_override(default_geometry='cylindrical', default_x_start=0.0d0, default_x_end=1.0d0)
    call initialise_grid()

    if (use_defaults) then
      cte_B02 = 0.1d0
      alpha = 2.0d0
      delta = 0.2d0

      k2 = 1.0d0
      k3 = 5.0d0

      resistivity = .true.
      use_fixed_resistivity = .true.
      fixed_eta_value = 1.0d-6

      viscosity = .true.
      viscosity_value = 1.0d-6
    end if

    B_field % B03 = 1.0d0

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      rho_field % rho0(i)       = 1.0d0 - (1.0d0 - delta) * r**2
      rho_field % d_rho0_dr(i)  = -2.0d0 * r * (1 - delta)

      if (int(alpha) == 0) then
        B_field % B02(i)          = cte_B02 * r / 2.0d0
        B_field % d_B02_dr(i)     = cte_B02 / 2.0d0
        p0(i)  = cte_B02**2 * (1 - r**2) / 4.0d0
        dp0(i) = -cte_B02**2 * r / 2.0d0
      else if (int(alpha) == 1) then
        B_field % B02(i)         = cte_B02 * r * (2.0d0 - r**2) / 4.0d0
        B_field % d_B02_dr(i)    = cte_B02 * (2.0d0 - 3.0d0 * r**2) / 4.0d0
        eta_field % dd_B02_dr(i) = -3.0d0 * cte_B02 * r / 2.0d0
        p0(i)  = cte_B02**2 * (5.0d0 / 12.0d0 - r**2 * (1.0d0 - 3.0d0 * r**2 / 4.0d0 + r**4 / 6.0d0)) / 4.0d0
        dp0(i) = - cte_B02**2 * r * (1.0d0 - 3.0d0 * r**2 / 2.0d0 + r**4 / 2.0d0) / 2.0d0
      else if (int(alpha) == 2) then
        B_field % B02(i)         = cte_B02 * r * (r**4 - 3.0d0 * r**2 + 3.0d0) / 6.0d0
        B_field % d_B02_dr(i)    = cte_B02 * (5.0d0 * r**4 - 9.0d0 * r**2 + 3.0d0) / 6.0d0
        eta_field % dd_B02_dr(i) = cte_B02 * r * (10.0d0 * r**2 - 9.0d0) / 3.0d0
        p0(i)  = cte_B02**2 * (47.0d0 / 360.0d0 - r**2 * (0.5d0 - 3.0d0 * r**2 / 4.0d0 &
                  + 5.0d0 * r**4 / 9.0d0 - 5.0d0 * r**6 / 24.0d0 + r**8 / 30.0d0)) / 2.0d0
        dp0(i) = -cte_B02**2 * r * (0.5d0 - 3.0d0 * r**2 / 2.0d0 + 5.0d0 * r**4 / 3.0d0 &
                  - 5.0d0 * r**6 / 6.0d0 + r**8 / 6.0d0)
      else
        call log_message("equilibrium: unexpected alpha value, try again with alpha = 0, 1, or 2", level='error')
      end if
      B_field % B0(i)      = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      T_field % T0(i)      = p0(i) / (rho_field % rho0(i))
      T_field % d_T0_dr(i) = dp0(i) / (rho_field % rho0(i)) &
                            - p0(i) * (rho_field % d_rho0_dr(i)) / (rho_field % rho0(i))**2
    end do

  end subroutine viscoresistive_tube_eq

end submodule smod_equil_viscoresistive_tube
