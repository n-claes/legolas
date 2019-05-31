module mod_setup_equilibrium
  use mod_global_variables
  implicit none

  real(dp), allocatable         ::  rho_0(:)
  real(dp), allocatable         ::  v_01(:), v_02(:), v_03(:)
  real(dp), allocatable         ::  T_0(:)
  real(dp), allocatable         ::  B_01(:), B_02(:), B_03(:)


contains

  subroutine initialise_equilibrium(grid, grid_gauss)

    real(dp), intent(in)  :: grid(gridpts)
    real(dp), intent(out) :: grid_gauss(4*gridpts)

    allocate(rho_0(4*gridpts))
    allocate(v_01(4*gridpts))
    allocate(v_02(4*gridpts))
    allocate(v_03(4*gridpts))
    allocate(T_0(4*gridpts))
    allocate(B_01(4*gridpts))
    allocate(B_02(4*gridpts))
    allocate(B_03(4*gridpts))

    call set_equilibrium(grid, grid_gauss)

  end subroutine initialise_equilibrium

  !> Sets the equilibrium on the nodes of the Gaussian quadrature
  subroutine set_equilibrium(grid, grid_gauss)

    real(dp), intent(in)  :: grid(gridpts)
    real(dp), intent(out) :: grid_gauss(4*gridpts)
    real(dp)              :: x_lo, x_hi, dx, xi(n_gauss)
    integer               :: i, j, idx

    ! Check for origin in cylindrical coordinates
    if (geometry == "cylindrical") then
      grid_gauss(1) = 1.0d-5
    end if

    ! Evaluate grid_gauss in nodes of Gaussian quadrature.
    ! An integral of f(x) in [a, b] can be approximated by
    ! 0.5*(b-a) * SUM[i from 1 -> n] ( wi * f( 0.5*(b-a)*xi + 0.5*(b-a)) )
    ! where wi and xi are the weights and nodes at i (so 1 to 4 here).
    ! Hence we need the gridpoints equal to the evaluation points of f(x).
    do i = 1, gridpts - 1
      x_lo = grid(i)
      x_hi = grid(i + 1)
      dx   = x_hi - x_lo

      do j = 1, n_gauss
        xi(j) = 0.5 * dx * gaussian_nodes(j) + 0.5*(x_lo + x_hi)
        idx   = (i - 1)*n_gauss + j
        grid_gauss(idx) = xi(j)
      end do
    end do

    ! Temporary initialisations
    rho_0 = 1.0d0
    v_01  = 1.0d0
    v_02  = 1.0d0
    v_03  = 1.0d0
    T_0   = 1.0d0
    B_01  = 1.0d0
    B_02  = 1.0d0
    B_03  = 1.0d0

  end subroutine set_equilibrium

  subroutine equilibrium_clean()
    deallocate(rho_0)
    deallocate(v_01)
    deallocate(v_02)
    deallocate(v_03)
    deallocate(T_0)
    deallocate(B_01)
    deallocate(B_02)
    deallocate(B_03)
  end subroutine equilibrium_clean


end module mod_setup_equilibrium
