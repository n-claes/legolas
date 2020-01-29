!
! SUBMODULE: smod_equil_gravito_acoustic
!
! DESCRIPTION:
!> Submodule defining an exponentially stratified hydrodynamic equilibrium
!! in Cartesian geometry with a constant gravity term included.
!! From Magnetohydrodynamics (2019), sec. 7.2.3
submodule (mod_equilibrium) smod_equil_gravito_acoustic
  implicit none

contains

  module subroutine gravito_acoustic_eq()
    real(dp)  :: x, g, rho0, p0, alpha
    integer   :: i

    geometry = 'Cartesian'
    x_start = 0.0d0
    x_end   = 1.0d0
    call initialise_grid()

    external_gravity = .true.

    k2 = dpi
    k3 = dpi

    !! Parameters
    g     = 0.5d0
    rho0  = 1.0d0
    p0    = 1.0d0
    alpha = rho0 * g / p0

    !! Equilibrium
    T_field % T0      = p0 / rho0
    grav_field % grav = g

    do i = 1, gauss_gridpts
      x = grid_gauss(i)

      !! Equilibrium
      rho_field % rho0(i) = rho0 * exp(-alpha*x)

      !! Derivative
      rho_field % d_rho0_dr(i) = -alpha * (rho_field % rho0(i))
    end do
    
  end subroutine gravito_acoustic_eq

end submodule smod_equil_gravito_acoustic
