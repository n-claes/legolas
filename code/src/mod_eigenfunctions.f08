module mod_eigenfunctions
  use mod_global_variables
  implicit none

  real(dp), allocatable         :: X_grid(:)
  complex(dp), allocatable      :: eigenfunctions(:, :)

contains

  subroutine initialise_eigenfunctions()
    allocate(X_grid(eigenf_gridpts))
    allocate(eigenfunctions(eigenf_gridpts, matrix_gridpts))

    X_grid = 0.0d0
    eigenfunctions = (0.0d0, 0.0d0)
  end subroutine initialise_eigenfunctions

  subroutine get_all_eigenfunctions(omegas, vr)
    use mod_io

    complex(dp), intent(in)     :: omegas(matrix_gridpts)
    complex(dp), intent(in)     :: vr(matrix_gridpts, matrix_gridpts)

    complex(dp)                 :: Y(eigenf_gridpts)
    integer                     :: i, j, var_idx

    character(len=3)            :: vars(8)

    vars = [character(3) :: 'rho', 'v1', 'v2', 'v3', 'T', 'a1', 'a2', 'a3']

    !! Calculate all eigenfunctions for rho
    do var_idx = 1, 8

      do i = 1, matrix_gridpts
        ! Every eigenfunction is stored in the column of 'eigenfunctions'.
        ! That way columns correspond with the eigenvalues and eigenvectors.
        ! So first column is first eigenfunction, second column is second, etc.
        call calculate_eigenfunction(trim(vars(var_idx)), vr(:, i), Y)
        eigenfunctions(:, i) = Y(:)
      end do
      
    end do
    stop

  end subroutine get_all_eigenfunctions

  !> Subroutine to calculate the eigenfunction of a given variable, taking
  !! the finite elements into account.
  !! @param[in] var   The variable (rho, v1, ...) for the eigenfunction
  !! @param[in] eigenvector The eigenvector for which to calculate the
  !!                        eigenfunction. This is one of the columns in
  !!                        'vr' from the solve_QR subroutine in mod_solvers
  !! @param[out] Y  The eigenfunction corresponding to the requested variable
  !!                and eigenvector
  subroutine calculate_eigenfunction(var, eigenvector, Y)
    use mod_grid

    character(len=*), intent(in)  :: var
    complex(dp), intent(in)       :: eigenvector(matrix_gridpts)
    complex(dp), intent(out)      :: Y(eigenf_gridpts)

    integer                       :: idx1, idx2, grid_idx
    real(dp)                      :: r, r_lo, r_hi
    real(dp)                      :: spline(4)
    character(len=3)              :: spline_type

    integer                       :: i

    !\TODO: re-transform variables


    if (var == 'rho') then
      idx1 = 1
      spline_type = "qua"
    else if (var == 'v1') then
      idx1 = 3
      spline_type = "cub"
    else if (var == 'v2') then
      idx1 = 5
      spline_type = "qua"
    else if (var == 'v3') then
      idx1 = 7
      spline_type = "qua"
    else if (var == 'T') then
      idx1 = 9
      spline_type = "qua"
    else if (var == 'a1') then
      idx1 = 11
      spline_type = "qua"
    else if (var == 'a2') then
      idx1 = 13
      spline_type = "cub"
    else if (var == 'a3') then
      idx1 = 15
      spline_type = "cub"
    else
      write(*, *) "Eigenfunction variable not known."
      write(*, *) "Currently set on: ", var
      stop
    end if

    idx2     = idx1 + 1
    grid_idx = 1

    !! Eigenfunctions are saved in each interval [grid(i), grid(i+1)] at
    !! the left edge, centre and right edge

    ! Contribution from first gridpoint, left edge
    r_lo = grid(1)
    r_hi = grid(2)
    r    = grid(1)

    X_grid(grid_idx) = r

    call get_spline(var, r, r_lo, r_hi, spline)
    Y(grid_idx) = Y(grid_idx) + eigenvector(idx1)*spline(2) &
                              + eigenvector(idx2)*spline(4)
    Y(grid_idx) = Y(grid_idx) + eigenvector(idx1 + dim_subblock)*spline(1) &
                              + eigenvector(idx2 + dim_subblock)*spline(3)

    ! Contribution from other gridpoints
    do i = 1, gridpts-1
      ! Contribution from centre
      r_lo = grid(i)
      r_hi = grid(i+1)
      r    = 0.5 * (r_lo + r_hi)
      ! save gridpoint
      grid_idx = grid_idx + 1
      X_grid(grid_idx) = r

      call get_spline(var, r, r_lo, r_hi, spline)
      Y(grid_idx) = Y(grid_idx) + eigenvector(idx1)*spline(2) &
                                + eigenvector(idx2)*spline(4)
      Y(grid_idx) = Y(grid_idx) + eigenvector(idx1 + dim_subblock)*spline(1) &
                                + eigenvector(idx2 + dim_subblock)*spline(3)

      ! Contribution from end point
      r = r_hi
      grid_idx = grid_idx + 1
      X_grid(grid_idx) = r

      call get_spline(var, r, r_lo, r_hi, spline)
      Y(grid_idx) = Y(grid_idx) + eigenvector(idx1)*spline(2) &
                                + eigenvector(idx2)*spline(4)
      Y(grid_idx) = Y(grid_idx) + eigenvector(idx1 + dim_subblock)*spline(1) &
                                + eigenvector(idx2 + dim_subblock)*spline(3)

      ! Increment indices of eigenvector elements
      idx1 = idx1 + dim_subblock
      idx2 = idx2 + dim_subblock
    end do

  end subroutine calculate_eigenfunction

  subroutine get_spline(var, r, r_lo, r_hi, spline)
    use mod_spline_functions

    character(len=*), intent(in)    :: var
    real(dp), intent(in)            :: r, r_lo, r_hi
    real(dp), intent(out)           :: spline(4)

    if ((var == "v1") .or. (var == "a2") .or. (var == "a3")) then
      call cubic_factors(r, r_lo, r_hi, spline)
    else
      call quadratic_factors(r, r_lo, r_hi, spline)
    end if

  end subroutine get_spline


  subroutine eigenfunctions_clean()
    deallocate(X_grid)
    deallocate(eigenfunctions)
  end subroutine eigenfunctions_clean


end module mod_eigenfunctions
