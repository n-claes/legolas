!
! SUBMODULE: smod_equil_coronal_flux_tube
!
! DESCRIPTION:
!> Submodule defining an equilibrium of a two-layer plasma.
!! Obtained from B. Roberts, MHD Waves in the Solar Atmosphere (2019).
!! Uses coronal conditions and should reproduce figure 6.7 (p. 173) in multiruns.
submodule(mod_equilibrium) smod_equil_coronal_flux_tube
  implicit none

contains

  module subroutine coronal_flux_tube_eq()
    use mod_global_variables, only: dp_LIMIT, gamma
    use mod_equilibrium_params, only: cte_rho0, cte_p0

    real(dp)  :: r, rho_e, p_e, B_0, B_e, r0
    integer   :: i

    call allow_geometry_override(default_geometry='cylindrical', default_x_start=0.0d0, default_x_end=1.0d0)
    call initialise_grid()

    if (use_defaults) then
      cte_rho0 = 1.0d0
      cte_p0 = 1.0d0
      r0 = 0.5d0

      k2 = 0.0d0
      k3 = 2.0d0
    end if

    rho_e = 4.0d0 * (2.0d0 * gamma + 1.0d0) * cte_rho0 / (50.0d0 * gamma + 1.0d0)
    p_e = (2.0d0 * gamma + 1.0d0) * cte_p0 / (50.0d0 * gamma + 1.0d0)
    B_0 = 2.0d0 * sqrt(gamma * cte_p0)
    B_e = 10.0d0 * sqrt(gamma * cte_p0 * (2.0d0 * gamma + 1.0d0) / (50.0d0 * gamma + 1.0d0))

    if (r0 > x_end) then
      error stop "equilibrium: inner cylinder radius r0 > x_end"
    else if (r0 < x_start) then
      error stop "equilibrium: inner cylinder radius r0 < x_start"
    end if

    ! check pressure balance
    if (abs(cte_p0 + 0.5d0 * B_0**2 - p_e - 0.5d0 * B_e**2) > dp_LIMIT) then
      error stop "equilibrium: total pressure balance not satisfied"
    end if

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      if (r > r0) then
        rho_field % rho0(i) = rho_e
        B_field % B02(i) = 0.0d0
        B_field % B03(i) = B_e
        T_field % T0(i) = p_e / rho_e
      else
        rho_field % rho0(i) = cte_rho0
        B_field % B02(i) = 0.0d0
        B_field % B03(i) = B_0
        T_field % T0(i) = cte_p0 / cte_rho0
      end if
    end do

  end subroutine coronal_flux_tube_eq

end submodule smod_equil_coronal_flux_tube
