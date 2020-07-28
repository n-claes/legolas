! =============================================================================
!> This submodule defines Rayleigh-Taylor instabilities in rotating theta pinches.
!! The straight cylinder approximation is used with a constant angular frequency.
!! Density and pressure profiles decrease over the domain, with a uni-directional
!! increasing magnetic field profile. Mode numbers \(k = 0\) correspond to
!! HD Rayleigh-Taylor instabilities, while \( k \neq 0 \) represent MHD RTIs.
!! The geometry is hardcoded to <tt>'cylindrical'</tt>, the domain is forced to
!! \(0 <= r <= 1\) through division by <tt>x_end</tt>.
!!
!! This equilibrium is taken from section IV in
!! _Goedbloed, J. P. "The Spectral Web of stationary plasma equilibria.
!! II. Internal modes." Physics of Plasmas 25.3 (2018): 032110_.
!! and also appears in section 13.4, figs. 13.11 to 13.15 in
!! _Goedbloed, H., Keppens, R., & Poedts, S. (2019). Magnetohydrodynamics of Laboratory
!! and Astrophysical Plasmas. Cambridge University Press._
!! [DOI](http://doi.org/10.1017/9781316403679).
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 1
!! - <tt>k3</tt> = 0 : so HD RTI
!! - <tt>cte_rho0</tt> = 1 : maximum density value.
!! - <tt>alpha</tt> = 2 : represents the stretching parameter.
!! - <tt>delta</tt> = 1/6 : represents the magnetic field deviation parameter.
!! - <tt>r0</tt> = 0 : represents the normalised radius \(x_0\) at maximum density.
!!
!! and can all be changed in the parfile. @endnote
! SUBMODULE: smod_equil_rotating_theta_pinch
!
! DESCRIPTION:
! Submodule defining Rayleigh-Taylor instabilities in a cylindrical geometry.
! Obtained from Goedbloed, Phys. Plasmas 25, 032110 (2018), Fig. 9, 11
! Also appears in Magnetohydrodynamics (2019), Fig. 13.12, 13.14
submodule (mod_equilibrium) smod_equil_RTI_theta_pinch
  implicit none

contains

  module subroutine RTI_theta_pinch_eq()
    use mod_equilibrium_params, only: cte_rho0, cte_p0, alpha, delta, r0

    real(dp)      :: B_inf, bigO, a
    real(dp)      :: r, x, fx, dfx
    integer       :: i

    flow = .true.
    geometry = 'cylindrical'
    call allow_geometry_override(default_x_start=0.0d0, default_x_end=1.0d0)
    call initialise_grid()

    if (use_defaults) then
      cte_rho0 = 1.0d0
      alpha = 2.0d0
      delta = 1.0d0 / 6.0d0
      r0 = 0.0d0

      k2 = 1.0d0
      k3 = 0.0d0
    end if

    a = x_end - x_start
    cte_p0 = 0.5d0 * (1.0d0 - delta)**2
    B_inf = a * sqrt(cte_rho0)
    bigO = alpha * sqrt(2.0d0 * delta * (1.0d0 - delta))

    do i = 1, gauss_gridpts
      r = grid_gauss(i)
      x = r / a
      fx = alpha**2 * (x**2 - r0**2)
      dfx = alpha**2 * 2.0d0 * x / a

      rho_field % rho0(i) = cte_rho0 / cosh(fx)**2
      v_field % v02(i) = bigO * r
      B_field % B03(i) = B_inf * (delta + (1.0d0 - delta) * tanh(fx))
      B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      T_field % T0(i) = cte_p0 / cte_rho0

      rho_field % d_rho0_dr(i) = -2.0d0 * cte_rho0 * dfx * tanh(fx) / cosh(fx)**2
      v_field % d_v02_dr(i) = bigO
      B_field % d_B03_dr(i) = B_inf * (1.0d0 - delta ) * dfx / cosh(fx)**2
    end do
  end subroutine RTI_theta_pinch_eq

end submodule smod_equil_RTI_theta_pinch
