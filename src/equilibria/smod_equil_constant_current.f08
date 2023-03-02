! =============================================================================
!> This submodule defines an equilibrium in cylindrical geometry with
!! a constant axial current. The geometry can be overridden using the parfile.
!!
!! This equilibrium is taken from
!! _Kerner, W. (1989). Large-scale complex eigenvalue problems.
!!  Journal of Computational Physics, 85(1), 1-85_.
!!
!! @note Default values are given by
!!
!! - <tt>k2</tt> = -2
!! - <tt>k3</tt> = 0.2
!! - <tt>j0</tt> = 0.125 : used to set the current.
!! - <tt>cte_rho0</tt> = 1 : used to set the density.
!! - <tt>cte_B03</tt> = 1 : used to set the Bz value.
!!
!! and can all be changed in the parfile. @endnote
submodule (mod_equilibrium) smod_equil_constant_current
  implicit none

contains

  !> Sets the equilibrium.
  module procedure constant_current_eq
    use mod_equilibrium_params, only: j0, cte_rho0, cte_B03

    real(dp) :: r
    real(dp), allocatable :: p_x(:), dp_x(:)
    integer :: i, gauss_gridpts

    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%grid%set_geometry("cylindrical")
      call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)
      k2 = -2.0d0
      k3 = 0.2d0

      j0 = 0.125d0
      cte_rho0 = 1.0d0
      cte_B03 = 1.0d0
    end if ! LCOV_EXCL_STOP
    call initialise_grid(settings)

    gauss_gridpts = settings%grid%get_gauss_gridpts()
    allocate(p_x(gauss_gridpts))
    allocate(dp_x(gauss_gridpts))

    rho_field % rho0  = cte_rho0
    B_field % B03     = cte_B03

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      B_field % B02(i)  = 0.5d0 * j0 * r
      B_field % B0(i)   = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      p_x(i)            = 0.25d0 * j0**2 * (1.0d0 - r**2)
      T_field % T0(i)   = p_x(i) / (rho_field % rho0(i))

      B_field % d_B02_dr(i) = 0.5d0 * j0
      dp_x(i)               = -0.5d0 * j0**2 * r
      T_field % d_T0_dr(i) = ( &
        dp_x(i) * (rho_field % rho0(i)) - (rho_field % d_rho0_dr(i)) * p_x(i) &
      ) / (rho_field % rho0(i))**2
    end do
    deallocate(p_x)
    deallocate(dp_x)
  end procedure constant_current_eq

end submodule smod_equil_constant_current
