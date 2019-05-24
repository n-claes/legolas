module mod_setup_equilibrium
  implicit none

  real, allocatable         ::  rho_0(:)
  real, allocatable         ::  v_0(:, :, :)
  real, allocatable         ::  T_0(:)
  real, allocatable         ::  B_0(:, :, :)


contains

  subroutine initialise_equilibrium(grid)
    use mod_global_variables

    double precision, intent(in)  ::  grid(gridpts)

    allocate(rho_0(gridpts))
    allocate(v_0(gridpts, gridpts, gridpts))
    allocate(T_0(gridpts))
    allocate(B_0(gridpts, gridpts, gridpts))

    call set_equilibrium(grid)

  end subroutine initialise_equilibrium

  !> Sets the equilibrium on the nodes of the Gaussian quadrature
  subroutine set_equilibrium(grid)
    use mod_global_variables

    double precision, intent(in)  :: grid(gridpts)
    double precision              :: grid_gauss(4*gridpts), xi(n_gauss)
    double precision              :: x_lo, x_hi, dx
    integer                       :: i, j, idx

    ! Fill odd indices of grid_gauss with the grid
    do i = 1, gridpts
      idx = (i - 1)*2 + 1
      grid_gauss(idx) = grid(i)
    end do
    ! Check for origin in cylindrical coordinates
    if (geometry == "cylindrical") then
      grid_gauss(1) = 1.0d-5
    end if

    ! Fill even indices of grid_gauss with midpoints between grid points
    do i = 1, gridpts-1
      idx = i*2
      grid_gauss(idx) = 0.5 * (grid(i) + grid(i+1))
    end do

    ! Evaluate grid in nodes of Gaussian quadrature
    do i = 1, gridpts - 1
      x_lo = grid(i)
      x_hi = grid(i + 1)
      dx   = x_hi - x_lo

      do j = 1, n_gauss
        xi(j) = x_lo + gaussian_nodes(j)*dx
        idx   = (i - 1)*n_gauss + j
        grid_gauss(idx) = xi(j)
      end do
    end do
    print*,grid
    print*,grid_gauss

    ! TODO: Why does ledaflowJD (line 850) use other node points?
    !       Maybe interval-related? [0, 1] instead of [-1, 1]?
    ! TODO: This method does not fill the last 4 points of grid_gauss?

    ! Temporary initialisations
    do i = 1, gridpts
      rho_0(i)     = 1.0d0
      v_0(i, i, i) = 1.0d0
      T_0(i)       = 1.0d0
      B_0(i, i, i) = 1.0d0
    end do

  end subroutine set_equilibrium

  subroutine equilibrium_clean()
    deallocate(rho_0)
    deallocate(v_0)
    deallocate(T_0)
    deallocate(B_0)
  end subroutine equilibrium_clean


end module mod_setup_equilibrium
