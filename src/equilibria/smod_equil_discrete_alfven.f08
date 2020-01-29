!
! SUBMODULE: smod_equil_discrete_alfven
!
! DESCRIPTION:
!> Submodule defining an equilibrium for non-adiabatic discrete Alfvén waves in
!! a cylindrical geometry.
!! Obtained from Keppens et al., Solar Physics 144, 267 (1993)
submodule (mod_equilibrium) smod_equil_discrete_alfven
  implicit none

contains

  module subroutine discrete_alfven_eq()
    real(dp)      :: r, x, j0, nu, nu_g, nu_l, d
    integer       :: i

    geometry = 'cylindrical'
    x_start = 0.0d0
    x_end   = 1.0d0
    call initialise_grid()

    flow = .true.

    ! Parameters
    nu    = 2.0d0
    nu_g  = nu + 1.0d0
    nu_l  = nu - 1.0d0
    j0    = 0.125d0
    d     = 0.2d0

    k2 = 1.0d0
    k3 = 0.05d0

    do i = 1, gauss_gridpts
      r = grid_gauss(i)
      x = r / x_end

      ! Equilibrium
      rho_field % rho0(i) = 1.0d0 - (1.0d0-d) * x**2
      v_field % v03(i)    = j0 * (1.0d0-x**2)**nu / (rho_field % rho0(i))
      B_field % B02(i)    = j0 * (1.0d0 - (1.0d0-x**2)**nu_g) / (2.0d0*x*nu_g)
      B_field % B03(i)    = 1.0d0
      B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      T_field % T0(i)     = 1.0d0 / (rho_field % rho0(i))

      ! Derivatives
      rho_field % d_rho0_dr(i) = -2.0d0 * (1.0d0-d) * x
      T_field % d_T0_dr(i)     = - (rho_field % d_rho0_dr(i)) / (rho_field % rho0(i))**2
      v_field % d_v03_dr(i)    = (-2.0d0*nu*j0*x*(1.0d0-r**2)**nu_l * (rho_field % rho0(i)) &
                                - (rho_field % d_rho0_dr(i)) * j0*(1.0d0-x**2)**nu) &
                                / (rho_field % rho0(i))**2
      B_field % d_B02_dr(i)    = j0 * (1.0d0-x**2)**nu - j0 * (1.0d0-(1.0d0-x**2)**nu_g) &
                                / (2.0d0*x**2*nu_g)

      ! TODO: eta?
      eta_field % dd_B02_dr(i) = (j0 / 2.0d0*nu_g) * ( - 4.0d0*nu_g*(1.0d0-x**2)**nu / x &
                      + (2.0d0*nu_g*(1-x**2)**nu - 4.0d0*nu*nu_g*x**2*(1.0d0-x**2)**nu_l) / x &
                      + 2.0d0*(1.0d0*(1.0d0-x**2)**nu_g) / x**3)
    end do

  end subroutine discrete_alfven_eq

end submodule smod_equil_discrete_alfven
