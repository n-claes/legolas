!
! SUBMODULE: smod_equil_gravito_mhd
!
! DESCRIPTION:
!> Submodule defining an exponentially stratified medium in Cartesian geometry
!! with a constant gravity term included.
!! From Magnetohydrodynamics (2019) by Goedbloed, Keppens and Poedts, sec. 7.3.3
submodule (mod_equilibrium) smod_equil_gravito_mhd
  implicit none

contains

  module subroutine gravito_mhd_eq()
    use mod_equilibrium_params, only: g, cte_rho0, cte_p0, alpha, beta

    real(dp)  :: x, B0
    integer   :: i

    geometry = 'Cartesian'
    x_start = 0.0d0
    x_end   = 1.0d0
    call initialise_grid()

    external_gravity = .true.

    if (use_defaults) then
      k2 = dpi
      k3 = dpi

      cte_p0 = 1.0d0
      g = 0.5d0
      alpha = 20.0d0
    end if

    B0 = 1.0d0
    beta  = 2.0d0*cte_p0 / B0**2
    cte_rho0 = (alpha / g) * (cte_p0 + 0.5d0 * B0**2)

    !! Equilibrium
    T_field % T0 = cte_p0 / cte_rho0
    grav_field % grav = g

    do i = 1, gauss_gridpts
      x = grid_gauss(i)

      !! Equilibrium
      rho_field % rho0(i) = cte_rho0 * exp(-alpha*x)
      B_field % B03(i)    = B0 * exp(-0.5d0 * alpha * x)
      B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)

      !! Derivatives
      rho_field % d_rho0_dr(i) = -alpha * (rho_field % rho0(i))
      B_field % d_B03_dr(i)    = -0.5d0 * alpha * (B_field % B03(i))
    end do

  end subroutine gravito_mhd_eq

end submodule smod_equil_gravito_mhd
