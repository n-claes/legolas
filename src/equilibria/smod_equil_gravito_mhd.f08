!
! SUBMODULE: smod_equil_gravito_mhd
!
! DESCRIPTION:
!> Submodule defining an exponentially stratified medium in Cartesian geometry
!! with a constant gravity term included.
!! From Magnetohydrodynamics (2019) by Goedbloed, Keppens and Poedts, sec. 7.3.3 (p258-261)
submodule (mod_equilibrium) smod_equil_gravito_mhd
  implicit none

contains

  module subroutine gravito_mhd_eq()
    use mod_equilibrium_params, only: g, cte_rho0, cte_p0, alpha, beta

    real(dp)  :: x, B0, tay_point
    integer   :: i
    logical, save :: use_taylor

    call allow_geometry_override(default_geometry='Cartesian', default_x_start=0.0d0, default_x_end=1.0d0)
    external_gravity = .true.

    ! For multiruns, this should be false. For single runs you should use the Taylor approximation and
    ! modified x_start and x_end values
    use_taylor = .false.

    if (use_defaults) then
      ! Using default values and Taylor is fine
      use_taylor = .true.
      k2 = dpi
      k3 = dpi

      cte_p0 = 0.5d0
      g = 0.5d0
      alpha = 20.0d0
    end if

    B0 = 1.0d0
    beta  = 2.0d0*cte_p0 / B0**2
    cte_rho0 = (alpha / g) * (cte_p0 + 0.5d0 * B0**2)

    if (use_taylor) then
      ! For Taylor approximation the interval should be small since the approximation is linear (2nd order).
      ! Therefore we divide x_end by the alpha value
      x_end   = x_end / alpha

      if ((x_end - x_start) > 0.1d0) then
        write(*, *) ">> WARNING: Using a 2nd order (linear) Taylor approximation for the exponent, but slab size > 0.1."
        write(*, *) "            The approximation may no longer be valid, and may result in false imaginary eigenvalues &
                                &due to numerical errors. Proceed with care."
        write(*, *) ""
      end if
    end if

    call initialise_grid()

    !! Equilibrium
    T_field % T0 = cte_p0 / cte_rho0
    grav_field % grav = g

    ! Taylor approximation, 2nd order (hence linear in x, but approximates fine in small [x_start, x_end] interval)
    if (use_taylor) then
      ! Centre of the Taylor approximation, taken to be in the middle of the grid.
      tay_point = (x_start + x_end) / 2.0d0

      do i = 1, gauss_gridpts
        x = grid_gauss(i)

        rho_field % rho0(i) = -alpha * cte_rho0 * (x - tay_point) * exp(-alpha * tay_point) &
                              + cte_rho0 * exp(-alpha * tay_point)
        rho_field % d_rho0_dr(i) = alpha**2 * cte_rho0 * (x - tay_point) * exp(-alpha * tay_point) &
                                   - alpha * cte_rho0 * exp(-alpha * tay_point)

        B_field % B03(i) = -0.5d0 * B0 * alpha * (x - tay_point) * exp(-0.5d0 * alpha * tay_point) &
                           + B0 * exp(-0.5d0 * alpha * tay_point)
        B_field % d_B03_dr(i) = 0.25d0 * B0 * alpha**2 * (x - tay_point) * exp(-0.5d0 * alpha * tay_point) &
                                - 0.5d0 * B0 * alpha * exp(-0.5d0 * alpha * tay_point)
        B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      end do

    ! Full exponential prescription for the equilibrium configuration
    else
      do i = 1, gauss_gridpts
        x = grid_gauss(i)

        rho_field % rho0(i) = cte_rho0 * exp(-alpha*x)
        B_field % B03(i) = B0 * exp(-0.5d0 * alpha * x)

        rho_field % d_rho0_dr(i) = -alpha * (rho_field % rho0(i))
        B_field % d_B03_dr(i) = -0.5d0 * alpha * (B_field % B03(i))
        B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      end do
    end if

  end subroutine gravito_mhd_eq

end submodule smod_equil_gravito_mhd
