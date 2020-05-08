!
! SUBMODULE: smod_equil_constant_current
!
! DESCRIPTION:
!> Initialises equilibrium with a constant axial current in cylindrical geometry
!! From Kerner, J. Comput. Phys. 85 (1989), Fig. 4.6
submodule (mod_equilibrium) smod_equil_constant_current
  implicit none

contains

  module subroutine constant_current_eq()
    use mod_equilibrium_params, only: j0

    real(dp)  :: r
    real(dp)  :: p_x(gauss_gridpts), dp_x(gauss_gridpts)
    integer   :: i

    call allow_geometry_override(default_geometry='cylindrical', default_x_start=0.0d0, default_x_end=1.0d0)
    call initialise_grid()

    if (use_defaults) then
      k2 = -2.0d0
      k3 = 0.2d0

      j0 = 0.125d0
    end if

    !! Equilibrium
    rho_field % rho0  = 1.0d0
    B_field % B03     = 1.0d0

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      !! Equilibrium
      B_field % B02(i)  = 0.5d0 * j0 * r
      B_field % B0(i)   = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      p_x(i)            = 0.25d0 * j0**2 * (1.0d0 - r**2)
      T_field % T0(i)   = p_x(i) / (rho_field % rho0(i))

      !! Derivatives
      B_field % d_B02_dr(i) = 0.5d0 * j0
      dp_x(i)               = -0.5d0 * j0**2 * r
      T_field % d_T0_dr(i) = ( dp_x(i) * (rho_field % rho0(i)) - (rho_field % d_rho0_dr(i)) * p_x(i) ) &
                              / (rho_field % rho0(i))**2
    end do

  end subroutine constant_current_eq

end submodule smod_equil_constant_current
