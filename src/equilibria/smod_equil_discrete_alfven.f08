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
    real(dp)      :: r, j0, nu, nu_g, d, c
    real(dp)      :: p_x(gauss_gridpts), dp_x(gauss_gridpts)
    integer       :: i

    geometry = 'cylindrical'
    x_start = 0.0d0
    x_end   = 1.0d0
    call initialise_grid()

    thermal_conduction = .true.

    ! Parameters
    nu    = 2.0d0
    nu_g  = nu + 1.0d0
    j0    = 0.125d0
    d     = 0.2d0
    c     = 47.0d0*j0**2 / 720.0d0    ! implies T = 0 at r = 1

    k2 = 1.0d0
    k3 = 0.05d0

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      ! Equilibrium
      rho_field % rho0(i) = 1.0d0 - (1.0d0-d) * r**2
      B_field % B02(i)    = j0 * (1.0d0 - (1.0d0-r**2)**nu_g) / (2.0d0*r*nu_g)
      B_field % B03(i)    = 1.0d0
      B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      ! This pressure only holds for n = 2
      ! Obtained from the equilibrium condition (integration constant c)
      p_x(i)              = c - j0**2 * (r**10/(60.0d0) - (5.0d0/48.0d0)*r**8 &
                            + (5.0d0/18.0d0)*r**6 - (3.0d0/8.0d0)*r**4 + r**2/4.0d0)
      T_field % T0(i)     = p_x(i) / (rho_field % rho0(i))

      ! Derivatives
      rho_field % d_rho0_dr(i) = 2.0d0 * (d-1.0d0) * r
      ! Again, only for nu = 2
      dp_x(i)                  = -j0**2 * (r**9/6.0d0 - (5.0d0/6.0d0)*r**7 &
                                  + (5.0d0/3.0d0)*r**5 - 1.5d0*r**3 + 0.5d0*r)
      T_field % d_T0_dr(i)     = (dp_x(i) * (rho_field % rho0(i)) &
                                  - (rho_field % d_rho0_dr(i)) * p_x(i)) &
                                  / (rho_field % rho0(i))**2
      B_field % d_B02_dr(i)    = j0 * (1.0d0-r**2)**nu - j0 * (1.0d0-(1.0d0-r**2)**nu_g) &
                                / (2.0d0*r**2*nu_g)
    end do

  end subroutine discrete_alfven_eq

end submodule smod_equil_discrete_alfven
