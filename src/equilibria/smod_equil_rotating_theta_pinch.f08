!
! SUBMODULE: smod_equil_rotating_theta_pinch
! 
! DESCRIPTION:
!> Submodule defining Rayleigh-Taylor instabilities in a cylindrical geometry.
!! Obtained from Goedbloed, Phys. Plasmas 25, 032110 (2018), Fig. 9, 11
!! Also appears in Magnetohydrodynamics (2019), Fig. 13.12, 13.14
submodule (mod_equilibrium) smod_equil_rotating_theta_pinch
  implicit none
  
contains 
  
  module subroutine rotating_theta_pinch_eq()
    real(dp)      :: p0, alpha, r, rho0, B_inf
    real(dp)      :: x, fx, dfx, ddfx, a, d, x0, bigO
    integer       :: i

    geometry = 'cylindrical'
    x_start = 0.0d0
    x_end   = 1.0d0
    call initialise_grid()

    flow = .true.

    !! Parameters
    rho0  = 1.0d0
    alpha = 2.0d0
    a     = x_end
    d     = 0.1667d0
    B_inf = 1.0d0
    p0    = 0.5d0 * (1.0d0-d)**2 * B_inf**2
    x0    = 0.0d0
    bigO  = alpha * sqrt(2.0d0*d*(1-d)) * B_inf / (a*sqrt(rho0))

    k2 = 1.0d0
    k3 = 0.0d0

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      x     = r / a
      fx    = alpha**2 * (x**2 - x0**2)
      dfx   = 2.0d0 * alpha**2 * x / a
      ddfx  = 2.0d0 * alpha**2 / a**2

      ! Equilibrium
      rho_field % rho0(i) = rho0 / cosh(fx)**2
      v_field % v02(i) = bigO * r
      B_field % B03(i) = B_inf * (d + (1.0d0-d) * tanh(fx))
      B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      T_field % T0(i) = p0 / rho0

      ! Derivatives
      rho_field % d_rho0_dr(i)  = -2.0d0*rho0 * dfx * tanh(fx) / cosh(fx)**2
      v_field % d_v02_dr = bigO
      B_field % d_B03_dr(i) = B_inf * (1-d) * dfx / cosh(fx)**2

      ! TODO: eta?
      eta_field % dd_B03_dr(i)  = B_inf * (1-d) * (ddfx - 2.0d0*dfx**2*tanh(fx)) / cosh(fx)**2
    end do
    
  end subroutine rotating_theta_pinch_eq
  
end submodule smod_equil_rotating_theta_pinch
