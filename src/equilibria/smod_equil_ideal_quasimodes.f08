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
    use mod_equilibrium_params, only: j0, p1, r0

    real(dp)      :: r, n
    real(dp)      :: p_r(gauss_gridpts), dp_r(gauss_gridpts)
    integer       :: i

    geometry = 'cylindrical'
    x_start = 0.0d0
    x_end   = 1.0d0
    call initialise_grid()

    resistivity = .true.
    use_fixed_resistivity = .true.
    fixed_eta_value = 5.0d-5

    if (use_defaults) then
      j0 = 0.5d0
      n = 1.0d0
      r0 = 10 * dpi
      k2 = 2.0d0
    else
      n = p1
    end if

    k3 = 2.0d0 * dpi * n / r0

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      rho_field % rho0(i) = 1.0d0
      B_field % B02(i) = j0 * r * (2.0d0 - r**2) / 4.0d0
      B_field % B03(i) = 1.0d0
      B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      p_r(i) = j0**2 * (2.0d0 * (x_end**6 - r**6) - 9.0d0 * (x_end**4 - r**4) + 12.0d0 * (x_end**2 - r**2)) / 48.0d0
      T_field % T0(i) = p_r(i) / rho_field % rho0(i)

      dp_r(i) = j0**2 * r * (3.0d0 * r**2 - r**4 - 2.0d0) / 4.0d0
      B_field % d_B02_dr(i) = j0 * (2.0d0 - 3.0d0 * r**2) / 4.0d0
      eta_field % dd_B02_dr(i) = -3.0d0 * j0 * r / 2.0d0
      T_field % d_T0_dr(i) = (dp_r(i) * (rho_field % rho0(i)) - (rho_field % d_rho0_dr(i)) * p_r(i)) &
                              / (rho_field % rho0(i))**2

    end do

  end subroutine ideal_quasimodes_eq

end submodule smod_equil_ideal_quasimodes
