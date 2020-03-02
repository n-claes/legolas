!
! SUBMODULE: smod_equil_resistive_tearing
!
! DESCRIPTION:
!> Submodule defining resistive tearing modes in Cartesian geometry without flow.
!! From Advanced Magnetohydrodynamics by Goedbloed, Keppens and Poedts, page 159.
submodule (mod_equilibrium) smod_equil_resistive_tearing
  implicit none

contains

  module subroutine resistive_tearing_modes_eq()
    use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value
    use mod_equilibrium_params, only: alpha, beta

    real(dp)              :: x
    integer               :: i

    geometry = 'Cartesian'
    x_start = -0.5d0
    x_end   = 0.5d0
    call initialise_grid()

    resistivity = .true.
    use_fixed_resistivity = .true.
    fixed_eta_value = 0.0001d0

    if (use_defaults) then
      k2 = 0.49d0
      k3 = 0.0d0

      alpha = 4.73884d0
      beta  = 0.15d0
    end if

    do i = 1, gauss_gridpts
      x = grid_gauss(i)

      ! Equilibrium
      rho_field % rho0(i) = 1.0d0
      B_field % B02(i)    = sin(alpha * x)
      B_field % B03(i)    = cos(alpha * x)
      B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      ! n=1, kB=1, mu0=1 for T0
      T_field % T0(i)     = beta * (B_field % B0(i))**2 / (2.0d0)

      ! Derivatives
      B_field % d_B02_dr(i) = alpha * cos(alpha * x)
      B_field % d_B03_dr(i) = -alpha * sin(alpha * x)
      ! No d_T0_dr needed, as B0**2 is independent of r

      eta_field % dd_B02_dr(i) = -alpha**2 * sin(alpha * x)
      eta_field % dd_B03_dr(i) = -alpha**2 * cos(alpha * x)
    end do
  end subroutine resistive_tearing_modes_eq

end submodule smod_equil_resistive_tearing
