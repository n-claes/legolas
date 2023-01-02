! =============================================================================
!> This submodule defines flow driven instabilities in a Cartesian geometry.
!! This equilibrium can not be called explicitly from the parfile, but rather acts
!! as a "parent setup" for the Rayleigh-Taylor and Kelvin-Helmholtz submodules which
!! use this specific kind of equilibrium but with different parameters.
!! This submodule is called within its implicit children.
!!
!! This equilibrium is taken from section 13.2, p. 486 in
!! _Goedbloed, H., Keppens, R., & Poedts, S. (2019). Magnetohydrodynamics of Laboratory
!!  and Astrophysical Plasmas. Cambridge University Press._ [DOI](http://doi.org/10.1017/9781316403679).
submodule (mod_equilibrium) smod_equil_flow_driven_instabilities
  implicit none

contains

  module procedure flow_driven_instabilities_eq
    use mod_equilibrium_params, only: g, delta, theta, p1, p2, p3, tau, &
                                      p4, alpha, cte_rho0, cte_p0

    real(dp)    :: v0, v1, v2, phi0, x
    real(dp)    :: v_prof, phi_prof, p_prof
    integer     :: i

    call settings%grid%set_geometry("Cartesian")
    call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)
    call initialise_grid(settings)
    call settings%physics%enable_flow()
    v0 = p1
    v1 = p2
    v2 = p3
    phi0 = p4

    do i = 1, settings%grid%get_gauss_gridpts()
      x = grid_gauss(i)

      v_prof = v0 + v1*(x - 0.5d0) + v2 * sin(tau * (x - 0.5d0))
      phi_prof = phi0 + alpha * (x - 0.5d0)
      p_prof = cte_p0 - (x - 0.5d0 * delta * x**2) * g

      rho_field % rho0(i) = cte_rho0 * (1.0d0 - delta*x)
      v_field % v02(i)    = sin(theta) * v_prof
      v_field % v03(i)    = cos(theta) * v_prof
      B_field % B02(i)    = sin(phi_prof)
      B_field % B03(i)    = cos(phi_prof)
      B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      T_field % T0(i)     = p_prof / (rho_field % rho0(i))
      grav_field % grav(i) = g

      rho_field % d_rho0_dr(i) = -cte_rho0 * delta
      v_field % d_v02_dr(i)    = sin(theta) * (v1 + v2 * tau * cos(tau * (x - 0.5d0)))
      v_field % d_v03_dr(i)    = cos(theta) * (v1 + v2 * tau * cos(tau * (x - 0.5d0)))
      B_field % d_B02_dr(i)    = cos(phi_prof) * alpha
      B_field % d_B03_dr(i)    = -sin(phi_prof) * alpha
      T_field % d_T0_dr(i)     = (-g * cte_rho0 * (1.0d0 - delta * x)**2 + cte_rho0 * delta * p_prof) &
                                  / (rho_field % rho0(i))**2
    end do
  end procedure flow_driven_instabilities_eq

end submodule smod_equil_flow_driven_instabilities
