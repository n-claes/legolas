module mod_equilibrium_derivatives
  use mod_global_variables
  implicit none

  !! Derivatives of the equilibrium variable arrays.
  !! d_xxx_dr     means d(xxx)/dr
  !! d_xxx_yyy_dr means d(xxx/yyy)/dr

  real(dp), allocatable       :: d_rho0_dr(:)
  real(dp), allocatable       :: d_rB02_dr(:)
  real(dp), allocatable       :: d_B03_r_dr(:)
  real(dp), allocatable       :: d_rv02_dr(:)
  real(dp), allocatable       :: d_v03_dr(:)
  real(dp), allocatable       :: d_B03_dr(:)
  real(dp), allocatable       :: d_T0_dr(:)

contains

  subroutine initialise_equilibrium_derivatives()
    allocate(d_rho0_dr(4*gridpts))
    allocate(d_rB02_dr(4*gridpts))
    allocate(d_B03_r_dr(4*gridpts))
    allocate(d_rv02_dr(4*gridpts))
    allocate(d_v03_dr(4*gridpts))
    allocate(d_B03_dr(4*gridpts))
    allocate(d_T0_dr(4*gridpts))

    d_rho0_dr  = 0.0d0
    d_rB02_dr  = 0.0d0
    d_B03_r_dr = 0.0d0
    d_rv02_dr  = 0.0d0
    d_v03_dr   = 0.0d0
    d_B03_dr   = 0.0d0
    d_T0_dr    = 0.0d0

    call calculate_equilibrium_derivatives()

  end subroutine initialise_equilibrium_derivatives

  subroutine calculate_equilibrium_derivatives()
    use mod_grid
    use mod_equilibrium

    call get_array_derivative(grid_gauss, rho0_eq, d_rho0_dr)
    call get_array_derivative(grid_gauss, (grid_gauss * B02_eq), d_rB02_dr)
    call get_array_derivative(grid_gauss, (B03_eq / grid_gauss), d_B03_r_dr)
    call get_array_derivative(grid_gauss, (grid_gauss * v02_eq), d_rv02_dr)
    call get_array_derivative(grid_gauss, v03_eq, d_v03_dr)
    call get_array_derivative(grid_gauss, B03_eq, d_B03_dr)
    call get_array_derivative(grid_gauss, T0_eq, d_T0_dr)

  end subroutine calculate_equilibrium_derivatives

  subroutine equilibrium_derivatives_clean()
    if (allocated(d_rho0_dr)) then
      deallocate(d_rho0_dr)
    end if
    if (allocated(d_rB02_dr)) then
      deallocate(d_rB02_dr)
    end if
    if (allocated(d_B03_r_dr)) then
      deallocate(d_B03_r_dr)
    end if
    if (allocated(d_rv02_dr)) then
      deallocate(d_rv02_dr)
    end if
    if (allocated(d_v03_dr)) then
      deallocate(d_v03_dr)
    end if
    if (allocated(d_B03_dr)) then
      deallocate(d_B03_dr)
    end if
    if (allocated(d_T0_dr)) then
      deallocate(d_T0_dr)
    end if
  end subroutine equilibrium_derivatives_clean

end module mod_equilibrium_derivatives
