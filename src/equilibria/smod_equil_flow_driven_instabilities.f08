!
! SUBMODULE: smod_equil_flow_driven_instabilities
!
! DESCRIPTION:
!> Submodule defining flow-driven instabilities in Cartesian geometry.
!! From Advanced Magnetohydrodynamics by Goedbloed, Keppens and Poedts, p. 17 (sec. 13.2).
submodule (mod_equilibrium) smod_equil_flow_driven_instabilities
  implicit none

contains

  module subroutine flow_driven_instabilities_eq()
    real(dp)    :: rho0, delta, theta, v0, v1, v2, tau, phi0, alpha, B0, x, k0, g, p0
    real(dp)    :: v_x(gauss_gridpts), phi_x(gauss_gridpts), p_x(gauss_gridpts)
    integer     :: i, big

    geometry = 'Cartesian'
    x_start = 0.0d0
    x_end   = 1.0d0
    call initialise_grid()

    g = 15.0d0

    flow = .true.
    external_gravity = .true.

    ! Parameters
    k0 = 1.0d0
    delta = -5.0d0
    phi0 = -0.35d0 * dpi
    alpha = 0.0d0
    theta = 0.35d0 * dpi
    v0 = 0.2d0
    v1 = 0.6d0
    v2 = 0.0d0
    tau = 0.0d0

    k2 = 0.0d0
    k3 = k0

    rho0 = 1.0d0
    p0   = real(huge(big))  ! arbitrarily imposed such that T > 0
    B0   = 1.0d0
    grav_field % grav = g

    do i = 1, gauss_gridpts
      x = grid_gauss(i)
      v_x(i)   = v0 + v1*(x - 0.5d0) + v2*sin(tau * (x - 0.5d0))
      phi_x(i) = phi0 + alpha*(x - 0.5d0)
      p_x(i)   = p0 - (x - 0.5d0 * delta * x**2)*g

      !! Equilibrium
      rho_field % rho0(i) = rho0 * (1.0d0 - delta*x)
      v_field % v02(i)    = sin(theta) * v_x(i)
      v_field % v03(i)    = cos(theta) * v_x(i)
      B_field % B02(i)    = B0 * sin(phi_x(i))
      B_field % B03(i)    = B0 * cos(phi_x(i))
      B_field % B0(i)     = B0
      T_field % T0(i)     = p_x(i) / (rho_field % rho0(i))

      !! Derivatives
      rho_field % d_rho0_dr(i) = -rho0 * delta
      v_field % d_v02_dr(i)    = sin(theta) * (v1 + v2*cos(tau * (x - 0.5d0)) * tau)
      v_field % d_v03_dr(i)    = cos(theta) * (v1 + v2*cos(tau * (x - 0.5d0)) * tau)
      B_field % d_B02_dr(i)    = B0 * cos(phi_x(i)) * alpha
      B_field % d_B03_dr(i)    = -B0 * sin(phi_x(i)) * alpha
      T_field % d_T0_dr(i)     = (-g * rho0*(1.0d0 - delta * x)**2 + rho0*delta*p_x(i)) &
                                  / (rho_field % rho0(i))**2
    end do

  end subroutine flow_driven_instabilities_eq

end submodule smod_equil_flow_driven_instabilities
