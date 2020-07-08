! =============================================================================
!> This submodule defines an exponentially stratified medium in
!! Cartesian geometry with a constant gravity term and magnetic shear.
!! The geometry can be overridden in the parfile.
!!
!! This equilibrium is taken from section 12.1.3 in
!! _Goedbloed, H., Keppens, R., & Poedts, S. (2019). Magnetohydrodynamics of Laboratory
!!  and Astrophysical Plasmas. Cambridge University Press._ [DOI](http://doi.org/10.1017/9781316403679).
!! @note Default values are given by
!!
!! - <tt>k2</tt> = \(\pi\)
!! - <tt>k3</tt> = \(\pi\)
!! - <tt>cte_p0</tt> = 0.25 : pressure, used to set the plasma beta.
!! - <tt>g</tt> = 0.5 : gravitational constant.
!! - <tt>lambda</tt> = 0 : magnetic shear value.
!! - <tt>alpha</tt> = 20 : constant to constrain the density value.
!!
!! and can all be changed in the parfile. @endnote
submodule (mod_equilibrium) smod_equil_interchange_modes
  implicit none

contains

  !> Sets the equilibrium.
  module subroutine interchange_modes_eq()
    use mod_equilibrium_params, only: g, cte_rho0, cte_p0, alpha, beta, lambda

    real(dp)  :: x, B0
    integer   :: i

    call allow_geometry_override(default_geometry='Cartesian', default_x_start=0.0d0, default_x_end=1.0d0)
    call initialise_grid()

    if (use_defaults) then
      external_gravity = .true.

      k2 = dpi
      k3 = dpi

      cte_p0 = 0.25d0
      g = 0.5d0
      lambda = 0.0d0
      alpha = 20.0d0
    end if

    B0 = 1.0d0
    beta = 2.0d0*cte_p0 / B0**2
    cte_rho0 = (alpha / g) * (cte_p0 + 0.5d0 * B0**2)

    T_field % T0      = cte_p0 / cte_rho0
    grav_field % grav = g

    do i = 1, gauss_gridpts
      x = grid_gauss(i)

      rho_field % rho0(i) = cte_rho0 * exp(-alpha*x)
      B_field % B02(i)    = B0 * exp(-0.5d0 * alpha * x) * sin(lambda*x)
      B_field % B03(i)    = B0 * exp(-0.5d0 * alpha * x) * cos(lambda*x)
      B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)

      rho_field % d_rho0_dr(i) = -alpha * (rho_field % rho0(i))
      B_field % d_B02_dr(i)    = -0.5d0 * alpha * (B_field % B02(i)) &
                                       + lambda * (B_field % B03(i))
      B_field % d_B03_dr(i)    = -0.5d0 * alpha * (B_field % B03(i)) &
                                       - lambda * (B_field % B02(i))
    end do
  end subroutine interchange_modes_eq

end submodule smod_equil_interchange_modes
