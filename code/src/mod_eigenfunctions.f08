module mod_eigenfunctions
  use mod_global_variables
  implicit none

  real(dp), allocatable         :: X_grid(:)


contains

  !> Initialises the eigenfunction module, allocates X_grid.
  subroutine initialise_eigenfunctions()
    allocate(X_grid(eigenf_gridpts))

    X_grid = 0.0d0
  end subroutine initialise_eigenfunctions

  !> Subroutine to calculate all eigenfunctions from the right eigenvectors.
  !! These are consequently saved to the proper files.
  !! @param[in] vr    Array of right eigenvectors, where every column
  !!                  has the eigenvector corresponding to the eigenvalue
  !!                  in 'omega' at the same column index
  subroutine get_all_eigenfunctions(vr)
    use mod_types
    use mod_io
    use mod_grid

    complex(dp), intent(in)     :: vr(matrix_gridpts, matrix_gridpts)

    complex(dp)                 :: Y(eigenf_gridpts)
    complex(dp)                 :: current_eigenvec(matrix_gridpts)
    integer                     :: i, j, k, idx
    character(len=3)            :: vars(8)
    real(dp)                    :: r_mid, r_end
    complex(dp)                 :: fact

    type (eigenf_type)          :: var_eigenf

    allocate(var_eigenf % eigenfunctions(eigenf_gridpts, matrix_gridpts))

    vars = [character(3) :: 'rho', 'v1', 'v2', 'v3', 'T', 'a1', 'a2', 'a3']

    !! Calculate all eigenfunctions for every variable
    do j = 1, 8

      ! Index of the current variable
      var_eigenf % index = j
      ! Integer for file writing
      var_eigenf % write_out = eigenvectors_out + j
      ! Name of the current variable
      var_eigenf % var = trim(vars(j))
      ! Filename of the current variable, used for saving
      var_eigenf % savename = trim(trim(vars(j)) // '_eigenfunction')

      ! do i = 1, matrix_gridpts
      !   ! Every eigenfunction is stored in the column of 'eigenfunctions'.
      !   ! That way columns correspond with the eigenvalues and eigenvectors.
      !   call calculate_eigenfunction(trim(vars(j)), vr(:, i), Y)
      !   var_eigenf % eigenfunctions(:, i) = Y(:)
      ! end do

      do i = 1, matrix_gridpts
        current_eigenvec(:) = vr(:, i)
        X_grid(1) = grid(1)
        idx = (j - 1)*2 + 2
        Y(1) = current_eigenvec(idx)

        do k = 2, gridpts
          r_mid = (grid(k) + grid(k-1)) / 2.0d0
          r_end = grid(k)
          X_grid(2*(k-1))   = r_mid
          X_grid(2*(k-1)+1) = r_end

          idx = (k - 1)*dim_subblock + (j-1)*2 + 1

          call get_transform_factor(trim(vars(j)), r_mid, fact)
          Y(2*(k-1)) = fact * current_eigenvec(idx)

          call get_transform_factor(trim(vars(j)), r_end, fact)
          Y(2*(k-1)+1) = fact * current_eigenvec(idx+1)
        end do
        var_eigenf % eigenfunctions(:, i) = Y(:)
      end do

      call save_eigenfunctions(var_eigenf)
    end do

    call save_eigenf_grid(X_grid)
    deallocate(var_eigenf % eigenfunctions)

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
    complex(dp)                   :: fact

    integer                       :: i

    Y = (0.0d0, 0.0d0)

    if (var == 'rho') then
      idx1 = 1
    else if (var == 'v1') then
      idx1 = 3
    else if (var == 'v2') then
      idx1 = 5
    else if (var == 'v3') then
      idx1 = 7
    else if (var == 'T') then
      idx1 = 9
    else if (var == 'a1') then
      idx1 = 11
    else if (var == 'a2') then
      idx1 = 13
    else if (var == 'a3') then
      idx1 = 15
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

    call get_transform_factor(var, r, fact)
    call get_spline(var, r, r_lo, r_hi, spline)
    Y(grid_idx) = Y(grid_idx) + fact*eigenvector(idx1)*spline(2) &
                              + fact*eigenvector(idx2)*spline(4)
    Y(grid_idx) = Y(grid_idx) + fact*eigenvector(idx1+dim_subblock)*spline(1) &
                              + fact*eigenvector(idx2+dim_subblock)*spline(3)

    ! Contribution from other gridpoints
    do i = 1, gridpts-1
      ! Contribution from centre
      r_lo = grid(i)
      r_hi = grid(i+1)
      r    = 0.5 * (r_lo + r_hi)
      ! save gridpoint
      grid_idx = grid_idx + 1
      X_grid(grid_idx) = r

      call get_transform_factor(var, r, fact)
      call get_spline(var, r, r_lo, r_hi, spline)
      Y(grid_idx) = Y(grid_idx) + fact*eigenvector(idx1)*spline(2) &
                                + fact*eigenvector(idx2)*spline(4)
      Y(grid_idx) = Y(grid_idx) + fact*eigenvector(idx1+dim_subblock)*spline(1) &
                                + fact*eigenvector(idx2+dim_subblock)*spline(3)

      ! Contribution from end point
      r = r_hi
      grid_idx = grid_idx + 1
      X_grid(grid_idx) = r

      call get_transform_factor(var, r, fact)
      call get_spline(var, r, r_lo, r_hi, spline)
      Y(grid_idx) = Y(grid_idx) + fact*eigenvector(idx1)*spline(2) &
                                + fact*eigenvector(idx2)*spline(4)
      Y(grid_idx) = Y(grid_idx) + fact*eigenvector(idx1 + dim_subblock)*spline(1) &
                                + fact*eigenvector(idx2 + dim_subblock)*spline(3)

      ! Increment indices of eigenvector elements
      idx1 = idx1 + dim_subblock
      idx2 = idx2 + dim_subblock
    end do

  end subroutine calculate_eigenfunction

  !> Calls the right finite element spline depending on the variable passed.
  !! @param[in] var   The variable ('rho', 'v1', etc.)
  !! @param[in] r     The current position in the grid
  !! @param[in] r_lo  Left edge of the current grid interval
  !! @param[in] r_hi  Right edge of the current grid interval
  !! @param[out] spline The finite element basis function, depending on 'var'
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


  !> Subroutine to retrieve the transformation factor. For example, rho was
  !! rescaled by eps*rho = RHO, such that in this case 'fact' will be 1/eps.
  !! @param[in] var   The variable ('rho', 'v1', etc.)
  !! @param[in] r     The current position in the grid
  !! @param[out] fact Transformation factor depending on the variable
  subroutine get_transform_factor(var, r, fact)
    character(len=*), intent(in) :: var
    real(dp), intent(in)         :: r
    complex(dp), intent(out)     :: fact

    real(dp)                     :: eps

    if ((geometry == 'Cartesian') .or. (geometry == 'cartesian')) then
      eps = 1.0d0
    else if (geometry == 'cylindrical') then
      eps = r
    else
      write(*, *) "Geometry not defined correctly"
      write(*, *) "Currently set on: ", geometry
      stop
    end if

    !! \note: ic = (0.0d0, 1.0d0) and ir = (1.0d0, 0.0d0)
    if (var == 'rho') then
      fact = (1.0d0 / eps) * ir
    else if (var == 'v1') then
      fact = 1.0d0 / (ic * eps)
    else if (var == 'v2') then
      fact = 1.0d0 * ir
    else if (var == 'v3') then
      fact = (1.0d0 / eps) * ir
    else if (var == 'T') then
      fact = (1.0d0 / eps) * ir
    else if (var == 'a1') then
      fact = 1.0d0 / ic
    else if (var == 'a2') then
      fact = (1.0d0 / eps) * ir
    else if (var == 'a3') then
      fact = 1.0d0 * ir
    else
      write(*, *) "Var not defined correctly"
      stop
    end if

  end subroutine get_transform_factor

  !> Cleaning routine, deallocates the variables in this module.
  subroutine eigenfunctions_clean()
    deallocate(X_grid)
  end subroutine eigenfunctions_clean


end module mod_eigenfunctions
