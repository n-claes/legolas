!
! SUBMODULE: smod_equil_suydam_cluster
!
! DESCRIPTION:
!> Submodule defining Suydam cluster modes in a cylindrical geometry.
!! Obtained from Bondeson et al., Phys. Fluids 30 (1987).
submodule (mod_equilibrium) smod_equil_suydam_cluster
  implicit none

contains

  module subroutine suydam_cluster_eq()
    use mod_equilibrium_params, only: cte_v03, cte_p0, p1, alpha

    real(dp)      :: r
    real(dp)      :: J0, J1, DJ0, DJ1
    real(dp)      :: P0_eq(gauss_gridpts)
    integer       :: i

    call allow_geometry_override(default_geometry='cylindrical', default_x_start=0.0d0, default_x_end=1.0d0)
    call initialise_grid()

    flow = .true.

    if (use_defaults) then
      cte_v03 = 0.14d0
      cte_p0 = 0.05d0
      p1 = 0.1d0
      alpha = 2.0d0

      k2 = 1.0d0
      k3 = -1.2d0
    end if

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      J0  = bessel_jn(0, alpha * r)
      J1  = bessel_jn(1, alpha * r)
      DJ0 = -alpha * J1
      DJ1 = alpha * (0.5d0 * J0 - 0.5d0 * bessel_jn(2, alpha * r))

      ! Equilibrium
      rho_field % rho0(i) = 1.0d0
      v_field % v02(i)    = 0.0d0
      v_field % v03(i)    = cte_v03 * (1.0d0 - r**2)
      B_field % B02(i)    = J1
      B_field % B03(i)    = sqrt(1.0d0 - p1) * J0
      B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      P0_eq(i)            = cte_p0 + 0.5d0 * p1 * J0**2
      T_field % T0(i)     = P0_eq(i) / (rho_field % rho0(i))

      ! Derivatives
      T_field % d_T0_dr(i)  = p1 * J0 * DJ0
      v_field % d_v03_dr(i) = -2.0d0 * cte_v03 * r
      B_field % d_B02_dr(i) = DJ1
      B_field % d_B03_dr(i) = -alpha * sqrt(1.0d0 - p1) * J1
    end do

  end subroutine suydam_cluster_eq

end submodule smod_equil_suydam_cluster
