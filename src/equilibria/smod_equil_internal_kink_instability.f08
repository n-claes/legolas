!
! SUBMODULE: smod_equil_internal_kink_instability
!
! DESCRIPTION:
!> Submodule defining current-driven internal kink instability in cylindrical geometry.
!! Obtained from Goedbloed, Phys. Plasmas 25, 032110 (2018), Fig. 3, 5.
submodule (mod_equilibrium) smod_equil_internal_kink_instability
  implicit none

contains

  module subroutine internal_kink_eq()
    real(dp)      :: v_z0, p0, alpha, r, rho0, x
    real(dp)      :: J0, J1, J2, J3, DJ0, DJ1, DDJ0, DDJ1
    integer       :: i

    geometry = 'cylindrical'
    x_start = 0.0d0
    x_end   = 1.0d0         ! this is parameter a in the paper
    call initialise_grid()

    flow = .true.

    !! Parameters
    rho0  = 1.0d0
    v_z0  = 1.0d0
    p0    = 3.0d0
    alpha = 5.0d0 / x_end

    k2 = 1.0d0
    k3 = 0.16d0 * alpha

    do i = 1, gauss_gridpts
      r = grid_gauss(i)
      x = r / x_end

      J0   = bessel_jn(0, alpha * x)
      J1   = bessel_jn(1, alpha * x)
      J2   = bessel_jn(2, alpha * x)
      J3   = bessel_jn(3, alpha * x)
      DJ0  = -alpha * J1
      DJ1  = alpha * (0.5d0 * J0 - 0.5d0 * J2)
      DDJ0 = -alpha * DJ1
      DDJ1 = -alpha**2 * (0.75d0 * J1 - 0.25d0 * J3)

      ! Equilibrium
      rho_field % rho0(i) = rho0 * (1.0d0-x**2)
      v_field % v03(i)    = v_z0 * (1.0d0-x**2)
      B_field % B02(i)    = J1
      B_field % B03(i)    = J0
      B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      T_field % T0(i)     = p0 / (rho_field % rho0(i))

      ! Derivatives
      rho_field % d_rho0_dr(i) = -2.0d0*rho0*x
      T_field % d_T0_dr(i)     = 2.0d0*x * p0 / (rho0 * (1.0d0-x**2)**2)
      v_field % d_v03_dr(i)    = -2.0d0*v_z0*x
      B_field % d_B02_dr(i)    = DJ1
      B_field % d_B03_dr(i)    = DJ0
    end do

  end subroutine internal_kink_eq

end submodule smod_equil_internal_kink_instability
