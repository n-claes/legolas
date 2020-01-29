!
! SUBMODULE: smod_equil_tokamak_cyl
!
! DESCRIPTION:
!> Initializes equilibrium for a tokamak-like case, with finite resistivity
!! in cylindrical geometry.
!! Obtained from Kerner, J. Comput. Phys. 85, 1-85 (1989), Eq. (4.43)
submodule (mod_equilibrium) smod_equil_tokamak_cyl
  implicit none

contains

  module subroutine tokamak_cyl_eq()
    use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value, dp_LIMIT

    real(dp)  :: j0, nu, r
    integer   :: i

    geometry = 'cylindrical'
    ! Override values from par file
    x_start = 0.0d0
    x_end   = 1.0d0
    call initialise_grid()

    flow = .true.
    resistivity = .true.
    use_fixed_resistivity = .true.
    fixed_eta_value = 1.0d-4

    k2 = 2.0d0
    k3 = 1.0d0

    !! Parameters
    j0  = 1.0d0
    nu  = 1.0d0

    !! Equilibrium
    rho_field % rho0 = 1.0d0
    B_field % B03    = 1.0d0
    T_field % T0     = 1.0d0

    if (abs(nu) > dp_LIMIT) then
      do i = 1, gauss_gridpts
        r = grid_gauss(i)

        v_field % v03(i)      = j0 * (1.0d0-r**2)**nu / (rho_field % rho0(i))
        v_field % d_v03_dr(i) = 2.0d0 * r * nu * j0 * (1.0d0-r**2)**(nu-1.0d0) &
                                  / (rho_field % rho0(i))
      end do
    else
      v_field % v03  = j0 / (rho_field % rho0(1))
    end if
  end subroutine tokamak_cyl_eq

end submodule smod_equil_tokamak_cyl
