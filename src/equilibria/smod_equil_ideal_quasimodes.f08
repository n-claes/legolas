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

    real(dp)      :: r, j0, n, L, c
    integer       :: i

    geometry = 'cylindrical'
    x_start = 0.0d0
    x_end   = 1.0d0
    call initialise_grid()

    resistivity = .true.
    use_fixed_resistivity = .true.
    fixed_eta_value = 5.0d-5

    !! Parameters
    j0  = 0.5d0
    n   = 1.0d0
    L   = 10*dpi

    c = 5.0d0/192.0d0     ! implies T = 0 at r = 1 for j0 = 0.5

    k2 = 2.0d0
    k3 = 2.0d0*dpi*n/L

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      ! Equilibrium
      rho_field % rho0(i) = 1.0d0
      B_field % B02(i)    = j0 * (1.0d0 - (1.0d0-r**2)**2) / (4.0d0*r)
      B_field % B03(i)    = 1.0d0
      B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      ! T obtained by solving the equilibrium condition for p0 (integration constant c)
      T_field % T0(i)     = c - 0.25d0*j0**2*r**2 * (1.0d0 - 0.75d0*r**2 + r**4/(6.0d0))

      ! Derivative
      B_field % d_B02_dr(i) = -0.25d0*j0 * (1.0d0 - (1.0d0-r**2)**2) / r**2 &
                                + 0.5d0*j0 * (1.0d0-r**2)
      T_field % d_T0_dr(i)  = -0.5d0*j0**2*r * (1.0d0 - 1.5d0*r**2 + 0.5d0*r**4)
    end do

  end subroutine ideal_quasimodes_eq

end submodule smod_equil_ideal_quasimodes
