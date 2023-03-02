! =============================================================================
!> This submodule defines magneto-rotational instabilities in an accretion disk.
!! Due to the special nature of this equilibrium <tt>x_start</tt> is hardcoded
!! to one and can not be overridden in the parfile, the same goes for the geometry
!! which is hardcoded to <tt>'cylindrical'</tt>. The outer edge can be chosen freely.
!! This equilibrium is chosen in such a way that the angular rotation is of order unity,
!! implying Keplerian rotation. The thin-disk approximation is valid with small magnetic
!! fields, but still large enough to yield magneto-rotational instabilities.
!! Gravity is assumed to go like \(g \thicksim 1/r^2\).
!!
!! This equilibrium is taken from section V in
!! _Goedbloed, J. P. "The Spectral Web of stationary plasma equilibria.
!! II. Internal modes." Physics of Plasmas 25.3 (2018): 032110_.
!! and also appears in section 13.5, fig. 13.7 in
!! _Goedbloed, H., Keppens, R., & Poedts, S. (2019). Magnetohydrodynamics of Laboratory
!! and Astrophysical Plasmas. Cambridge University Press._
!! [DOI](http://doi.org/10.1017/9781316403679).
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 0
!! - <tt>k3</tt> = 70
!! - <tt>beta</tt> = 100 : parameter \(\beta = 2p_1/B_1^2\).
!! - <tt>tau</tt> = 1 : represents parameter \(\tau = \mu_1 = B_{\theta 1}/B_{z1}\)
!! - <tt>nu</tt> = 0.1 : represents parameter \(\nu = \epsilon = \sqrt{p_1}\)
!! - <tt>x_end</tt> = 2 : fixes the parameter \(\delta = x_{end}/x_{start}\)
!!
!! and can all be changed in the parfile. @endnote
submodule (mod_equilibrium) smod_equil_MRI
  implicit none

contains

  module procedure MRI_accretion_eq
    use mod_equilibrium_params, only: beta, tau, nu

    real(dp) :: r, Bth1, Bz1, p1, delta, epsilon, mu1, vth1
    real(dp) :: x_start, x_end
    real(dp), allocatable :: p0(:), dp_dr(:)
    integer :: i, gauss_gridpts

    call settings%grid%set_geometry("cylindrical")

    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%grid%set_grid_boundaries(1.0_dp, 2.0_dp)

      call settings%physics%enable_flow()
      call settings%physics%enable_gravity()
      k2 = 0.0d0
      k3 = 70.0d0
      beta = 100.0d0
      tau = 1.0d0
      nu = 0.1d0
    end if ! LCOV_EXCL_STOP
    call initialise_grid(settings)

    gauss_gridpts = settings%grid%get_gauss_gridpts()
    allocate(p0(gauss_gridpts))
    allocate(dp_dr(gauss_gridpts))
    x_start = settings%grid%get_grid_start()
    x_end = settings%grid%get_grid_end()

    mu1 = tau
    epsilon = nu

    delta = x_end / x_start
    p1 = epsilon**2
    Bz1 = sqrt(2.0d0 * p1 / (beta * (1.0d0 + mu1**2)))
    Bth1 = mu1 * Bz1
    vth1 = sqrt(1.0d0 - 2.5d0 * p1 - 0.25d0 * Bth1**2 - 1.25d0 * Bz1**2)

    do i = 1, settings%grid%get_gauss_gridpts()
      r = grid_gauss(i)

      p0(i) = p1 * r**(-2.5d0)
      dp_dr(i) = -2.5d0 * p1 * r**(-3.5d0)

      rho_field % rho0(i) = r**(-1.5d0)
      T_field % T0(i) = p0(i) / (rho_field % rho0(i))
      v_field % v02(i) = vth1 / sqrt(r)
      B_field % B02(i) = Bth1 * r**(-1.25d0)
      B_field % B03(i) = Bz1 * r**(-1.25d0)
      B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      grav_field % grav(i) = 1.0d0 / r**2

      rho_field % d_rho0_dr(i) = -1.5d0 * r**(-2.5d0)
      T_field % d_T0_dr(i) = ( &
        dp_dr(i) * (rho_field % rho0(i)) - (rho_field % d_rho0_dr(i)) * p0(i) &
      ) / (rho_field % rho0(i))**2
      v_field % d_v02_dr(i) = -0.5d0 * vth1 * r**(-1.5d0)
      B_field % d_B02_dr(i) = -1.25d0 * Bth1 * r**(-2.25d0)
      B_field % d_B03_dr(i) = -1.25d0 * Bz1 * r**(-2.25d0)
    end do
    deallocate(p0)
    deallocate(dp_dr)
  end procedure MRI_accretion_eq

end submodule smod_equil_MRI
