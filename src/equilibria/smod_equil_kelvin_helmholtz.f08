!
! SUBMODULE: smod_equil_kelvin_helmholtz
!
! DESCRIPTION:
! Submodule defining a Kelvin-Helmholtz instability in Cartesian geometry.
! From Miura et al., J. Geophys. Res. 87 (1982)
submodule (mod_equilibrium) smod_equil_kelvin_helmholtz
  implicit none

contains

  module subroutine kh_instability_eq()
    use mod_equilibrium_params, only: cte_p0, cte_v02, cte_v03, p1

    real(dp)    :: a, x
    integer     :: i

    geometry = 'Cartesian'
    ! if on interval [0,1], change velocity profile x --> x - 0.5
    x_start = -0.5d0
    x_end   = 0.5d0
    call initialise_grid()

    flow = .true.

    if (use_defaults) then
      k2 = 10.0d0
      k3 = 0.0d0

      a = 0.05d0
      cte_p0  = 3.6d0
      cte_v02 = 1.67d0
      cte_v03 = 0.0d0   ! so B is perpendicular to v
    else
      a = p1
    end if

    !! Filling equilibria
    rho_field % rho0 = 1.0d0
    B_field % B02    = 0.0d0
    B_field % B03    = 1.0d0
    B_field % B0     = sqrt((B_field % B02)**2 + (B_field % B03**2))

    do i = 1, gauss_gridpts
      x = grid_gauss(i)
      v_field % v02(i) = -cte_v02 * tanh(x / a)
      v_field % v03(i) = -cte_v03 * tanh(x / a)
      T_field % T0(i)  = cte_p0 / (rho_field % rho0(i))

      ! Derivatives
      v_field % d_v02_dr(i) = -cte_v02 / (a * cosh(x / a)**2)
      v_field % d_v03_dr(i) = -cte_v03 / (a * cosh(x / a)**2)
    end do

  end subroutine kh_instability_eq

end submodule smod_equil_kelvin_helmholtz
