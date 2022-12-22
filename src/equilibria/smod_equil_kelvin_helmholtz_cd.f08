! =============================================================================
!> This submodule defines an unperturbed magnetised jet model in cylindrical
!! geometry, giving rise to Kelvin-Helmholtz and current-driven instabilities.
!! The geometry is fixed for this problem; the cylinder wall is dependent
!! on the equilibrium parameters and is given by <tt>2rj</tt>.
!!
!! This equilibrium is taken from
!! _Baty, H., & Keppens, R. (2002). Interplay between Kelvin-Helmholtz and
!!  current-driven instabilities in jets. The Astrophysical Journal, 580(2), 800_.
!! @note Default values are given by
!!
!! - <tt>k2</tt> = -1
!! - <tt>V</tt> = 1.63 : amplitude of the velocity shear
!! - <tt>cte_p0</tt> = 1 : used to set the pressure.
!! - <tt>cte_rho0</tt> = 1 : used to set the density.
!! - <tt>Bz0</tt> = 0.25 : used to set Bz.
!! - <tt>rc</tt> = 0.5 : length for radial variation.
!! - <tt>rj</tt> = 1 : jet radius
!!
!! and can all be changed in the parfile. @endnote
!! @note The default setup is _HEL2_ in the original paper.
!!       For _HEL1_ you can set <tt>rc = 2</tt>. @endnote
submodule (mod_equilibrium) smod_equil_kelvin_helmholtz_cd
  implicit none

contains

  !> Sets the equilibrium.
  module procedure kh_cd_instability_eq
    use mod_equilibrium_params, only: V, cte_rho0, cte_p0, Bz0, rc, Bth0, rj

    real(dp)    :: r, a
    integer     :: i

    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%physics%enable_flow()

      V = 1.63d0
      cte_rho0 = 1.0d0
      cte_p0 = 1.0d0
      Bz0 = 0.25d0
      rj = 1.0d0

      rc    = 0.5d0
      k2  = -1.0d0
    end if ! LCOV_EXCL_STOP

    Bth0  = 0.4d0 * (rc**2+rj**2) / (rj*rc)
    a = 0.1d0 * rj
    k3  = dpi / rj

    call settings%grid%set_geometry("cylindrical")
    call settings%grid%set_grid_boundaries(0.0_dp, 2.0_dp * rj)
    call initialise_grid(settings)

    do i = 1, settings%grid%get_gauss_gridpts()
      r = grid_gauss(i)

      rho_field % rho0(i) = cte_rho0
      v_field % v03(i)    = (V/2.0d0) * tanh((rj-r)/a)
      B_field % B02(i)    = Bth0 * r*rc / (rc**2 + r**2 )
      B_field % B03(i)    = Bz0
      B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      T_field % T0(i)     = cte_p0 / (rho_field % rho0(i)) &
        - (Bth0**2/(2.0d0*(rho_field % rho0(i)))) * (1 - rc**4/(rc**2+r**2)**2)

      B_field % d_B02_dr(i) = Bth0 * rc * (rc**2-r**2) / (r**2+rc**2)**2
      v_field % d_v03_dr(i) = - (V/(2.0d0*a)) / cosh((rj-r)/a)**2
      T_field % d_T0_dr(i)  = - (2.0d0*Bth0**2/(rho_field % rho0(i))) &
        * rc**4*r / (r**2+rc**2)**3
    end do
  end procedure kh_cd_instability_eq

end submodule smod_equil_kelvin_helmholtz_cd
