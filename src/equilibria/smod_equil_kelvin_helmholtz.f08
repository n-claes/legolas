!
! SUBMODULE: smod_equil_kelvin_helmholtz
! 
! DESCRIPTION:
!> Submodule defining a Kelvin-Helmholtz instability in Cartesian geometry.
!! From Miura et al., J. Geophys. Res. 87 (1982)
submodule (mod_equilibrium) smod_equil_kelvin_helmholtz
  implicit none
  
contains
  
  module subroutine kh_instability_eq()
    real(dp)    :: a, x, p0, v0y, v0z
    integer     :: i

    geometry = 'Cartesian'
    ! if on interval [0,1], change velocity profile x --> x - 0.5
    x_start = -0.5d0
    x_end   = 0.5d0
    call initialise_grid()

    flow = .true.

    !! Parameters
    a   = 0.05d0
    P0  = 3.6d0
    v0y = 1.67d0
    v0z = 0.0d0   ! so B is perpendicular to v

    !! Filling equilibria
    rho_field % rho0 = 1.0d0
    B_field % B02 = 0.0d0
    B_field % B03 = 1.0d0
    B_field % B0 = sqrt((B_field % B02)**2 + (B_field % B03**2))

    k2 = 10.0d0
    k3 = 0.0d0

    do i = 1, gauss_gridpts
      x = grid_gauss(i)
      v_field % v02(i) = -v0y * tanh(x / a)
      v_field % v03(i) = -v0z * tanh(x / a)
      T_field % T0(i) = P0 / (rho_field % rho0(i))

      ! Derivatives
      v_field % d_v02_dr(i) = -v0y / (a * cosh(x / a)**2)
      v_field % d_v03_dr(i) = -v0z / (a * cosh(x / a)**2)
    end do
  
  end subroutine kh_instability_eq
  
end submodule smod_equil_kelvin_helmholtz
    
