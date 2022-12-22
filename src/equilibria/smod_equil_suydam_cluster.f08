! =============================================================================
!> This submodule defines a cylindrical equilibrium in which a Suydam surface
!! is present, such that this gives rises to Suydam cluster modes.
!! The geometry can be overriden using the parfile.
!!
!! This equilibrium is taken from
!! _Nijboer, R. J., Holst, B., Poedts, S., & Goedbloed, J. P. (1997).
!!  Calculating magnetohydrodynamic flow spectra. Computer physics communications, 106(1-2), 39-52_.
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 1
!! - <tt>k3</tt> = -1.2
!! - <tt>cte_rho0</tt> = 1 : used to set the density value.
!! - <tt>cte_p0</tt> = 0.05 : used to set the background pressure.
!! - <tt>cte_v02</tt> = 0 : sets a constant flow \(v_\theta\).
!! - <tt>cte_v03</tt> = 0.14 : prefactor for flow profile \(v_z\).
!! - <tt>p1</tt> = 0.1 : constant that appears in magnetic field and pressure components.
!! - <tt>alpha</tt> = 2 : used in the Bessel functions.
!!
!! and can all be changed in the parfile. @endnote
submodule (mod_equilibrium) smod_equil_suydam_cluster
  implicit none

contains

  !> Sets the equilibrium
  module procedure suydam_cluster_eq
    use mod_equilibrium_params, only: cte_rho0, cte_v02, cte_v03, cte_p0, p1, alpha

    real(dp) :: r
    real(dp) :: J0, J1, DJ0, DJ1
    real(dp), allocatable:: P0_eq(:)
    integer :: i, gauss_gridpts

    if (settings%equilibrium%use_defaults) then  ! LCOV_EXCL_START
      call settings%grid%set_geometry("cylindrical")
      call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)
      call settings%physics%enable_flow()

      cte_rho0 = 1.0d0
      cte_v02 = 0.0d0
      cte_v03 = 0.14d0
      cte_p0 = 0.05d0
      p1 = 0.1d0
      alpha = 2.0d0

      k2 = 1.0d0
      k3 = -1.2d0
    end if ! LCOV_EXCL_STOP
    call initialise_grid(settings)

    gauss_gridpts = settings%grid%get_gauss_gridpts()
    allocate(P0_eq(gauss_gridpts))
    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      J0  = bessel_jn(0, alpha * r)
      J1  = bessel_jn(1, alpha * r)
      DJ0 = -alpha * J1
      DJ1 = alpha * (0.5d0 * J0 - 0.5d0 * bessel_jn(2, alpha * r))

      rho_field % rho0(i) = cte_rho0
      v_field % v02(i)    = cte_v02
      v_field % v03(i)    = cte_v03 * (1.0d0 - r**2)
      B_field % B02(i)    = J1
      B_field % B03(i)    = sqrt(1.0d0 - p1) * J0
      B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      P0_eq(i)            = cte_p0 + 0.5d0 * p1 * J0**2
      T_field % T0(i)     = P0_eq(i) / (rho_field % rho0(i))

      T_field % d_T0_dr(i)  = p1 * J0 * DJ0
      v_field % d_v03_dr(i) = -2.0d0 * cte_v03 * r
      B_field % d_B02_dr(i) = DJ1
      B_field % d_B03_dr(i) = -alpha * sqrt(1.0d0 - p1) * J1
    end do
    deallocate(P0_eq)
  end procedure suydam_cluster_eq

end submodule smod_equil_suydam_cluster
