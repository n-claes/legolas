!> Submodule for user-defined equilibria.
!! Look at the examples in the equilibria subdirectory or consult the website for
!! more information.
submodule (mod_equilibrium) smod_user_defined
  implicit none

contains

  ! LCOV_EXCL_START <exclude this file from code coverage>
  module subroutine user_defined_eq()
    real(dp)    :: x
    integer     :: i

    ! either allow for an override through the parfile or hardcode it
    call allow_geometry_override( &
      default_geometry="Cartesian", default_x_start=0.0d0, default_x_end=1.0d0 &
    )
    call initialise_grid()

    k2 = 0.0d0
    k3 = 1.0d0

    ! additional physics
    flow = .true.

    ! set up the grid
    do i = 1, gauss_gridpts
      x = grid_gauss(i)

      ! Note: values that are not set/referenced are automatically set to zero.

      ! standard equilibrium quantities
      rho_field%rho0(i) = 1.0d0
      v_field%v03(i) = 1.0d0
      T_field%T0(i) = 1.0d0
      B_field%B03(i) = 1.0d0
      B_field%B0(i) = sqrt(B_field%B02(i)**2 + B_field%B03(i)**2)

      ! first order derivatives
      v_field%d_v03_dr(i) = 0.0d0

      ! in case of resistivity, set second order B derivatives
      eta_field%dd_B03_dr = 0.0d0
    end do

  end subroutine user_defined_eq
  ! LCOV_EXCL_STOP

end submodule smod_user_defined
