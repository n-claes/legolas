!
! SUBMODULE: smod_equil_ideal_quasimodes
!
! DESCRIPTION:
!> Submodule defining ideal quasimodes in a cylindrical geometry.
!! Obtained from Poedts & Kerner, Physical Review Letters 66 (22), (1991)
submodule (mod_equilibrium) smod_equil_ideal_quasimodes
  implicit none

contains

  module subroutine ideal_quasimodes_eq()
    use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value

    real(dp)      :: r, x, j0, n, L, beta
    integer       :: i

    geometry = 'cylindrical'
    x_start = 0.0d0
    x_end   = 1.6d0
    call initialise_grid()

    flow = .true.
    use_fixed_resistivity = .true.
    fixed_eta_value = 5.0d-5

    !! Parameters
    j0  = 0.5d0
    n   = 1.0d0
    L   = 10*dpi

    beta = 2.0d0

    k2 = 2.0d0
    k3 = 2.0d0*dpi*n/L

    do i = 1, gauss_gridpts
      r = grid_gauss(i)
      x = r / x_end

      ! Equilibrium
      rho_field % rho0(i) = 1.0d0
      v_field % v03(i)    = j0 * (1.0d0-x**2) / (rho_field % rho0(i))
      B_field % B03(i)    = 1.0d0
      B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      ! n=1, kB=1, mu0=1 for T0
      T_field % T0(i)     = beta * (B_field % B0(i))**2 / (2.0d0*(rho_field % rho0(i)))

      ! Derivative
      v_field % d_v03_dr(i)  = -2.0d0*j0*x / (x_end*(rho_field % rho0(i)))
    end do

  end subroutine ideal_quasimodes_eq

end submodule smod_equil_ideal_quasimodes
