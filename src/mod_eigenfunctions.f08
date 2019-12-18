module mod_eigenfunctions
  use mod_global_variables
  implicit none

  real(dp), allocatable         :: ef_grid(:)


contains

  !> Initialises the eigenfunction module, allocates X_grid.
  subroutine initialise_eigenfunctions()
    allocate(ef_grid(ef_gridpts))

    ef_grid = 0.0d0
  end subroutine initialise_eigenfunctions

  !> Subroutine to calculate all eigenfunctions from the right eigenvectors.
  !! These are consequently saved to the proper files.
  !! @param[in] vr    Array of right eigenvectors, where every column
  !!                  has the eigenvector corresponding to the eigenvalue
  !!                  in 'omega' at the same column index
  subroutine calculate_eigenfunctions(vr)
    use mod_types, only: ef_type
    use mod_output, only: eigenfunctions_tofile, ef_grid_tofile

    complex(dp), intent(in)     :: vr(matrix_gridpts, matrix_gridpts)

    complex(dp)                 :: Y(ef_gridpts)
    integer                     :: i, j
    character(len=3)            :: vars(8)
    type (ef_type)              :: var_ef

    allocate(var_ef % eigenfunctions(ef_gridpts, matrix_gridpts))

    vars = [character(3) :: 'rho', 'v1', 'v2', 'v3', 'T', 'a1', 'a2', 'a3']

    !! Calculate all eigenfunctions for every variable
    do j = 1, 8
      ! Index of the current variable
      var_ef % index = j
      ! Name of the current variable
      var_ef % var = trim(vars(j))
      ! Filename used for saving (eg. rho_eigenfunction)
      var_ef % savename = trim(trim(vars(j)) // '_eigenfunction')

      do i = 1, matrix_gridpts
        ! The eigenfunction for each omega is stored in the column of
        ! 'eigenfunctions' with the same index as omega.
        ! Column indices hence correspond to the eigenvalue at that index.
        call get_eigenfunction(trim(vars(j)), vr(:, i), Y)
        var_ef % eigenfunctions(:, i) = Y(:)
      end do

      call eigenfunctions_tofile(var_ef)
    end do

    call ef_grid_tofile(ef_grid, trim('ef_' // savename_grid))
    deallocate(var_ef % eigenfunctions)

  end subroutine calculate_eigenfunctions

  !> Subroutine to calculate the eigenfunction of a given variable, taking
  !! the finite elements into account.
  !! @param[in] var   The variable (rho, v1, ...) for the eigenfunction
  !! @param[in] eigenvector The eigenvector for which to calculate the
  !!                        eigenfunction. This is one of the columns in 'vr'
  !! @param[out] Y  The eigenfunction corresponding to the requested variable
  !!                and eigenvector
  subroutine get_eigenfunction(var, eigenvector, Y)
    use mod_grid, only: grid

    character(len=*), intent(in)  :: var
    complex(dp), intent(in)       :: eigenvector(matrix_gridpts)
    complex(dp), intent(out)      :: Y(ef_gridpts)

    integer                       :: idx1, idx2, grid_idx, i
    real(dp)                      :: r, r_lo, r_hi
    real(dp)                      :: spline(4)
    complex(dp)                   :: fact

    Y = (0.0d0, 0.0d0)

    select case(var)
      case('rho')
        idx1 = 1
      case('v1')
        idx1 = 3
      case('v2')
        idx1 = 5
      case('v3')
        idx1 = 7
      case('T')
        idx1 = 9
      case('a1')
        idx1 = 11
      case('a2')
        idx1 = 13
      case('a3')
        idx1 = 15
      case default
        write(*, *) "Eigenfunction variable not known."
        write(*, *) "Currently set on: ", var
        stop
    end select

    idx2     = idx1 + 1
    grid_idx = 1

    !! Eigenfunctions are saved in each interval [grid(i), grid(i+1)] at
    !! the left edge, centre and right edge

    ! Contribution from first gridpoint, left edge
    r_lo = grid(1)
    r_hi = grid(2)
    r    = grid(1)

    ef_grid(grid_idx) = r

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
      ef_grid(grid_idx) = r

      call get_transform_factor(var, r, fact)
      call get_spline(var, r, r_lo, r_hi, spline)
      Y(grid_idx) = Y(grid_idx) + fact*eigenvector(idx1)*spline(2) &
                                + fact*eigenvector(idx2)*spline(4)
      Y(grid_idx) = Y(grid_idx) + fact*eigenvector(idx1+dim_subblock)*spline(1) &
                                + fact*eigenvector(idx2+dim_subblock)*spline(3)

      ! Contribution from end point
      r = r_hi
      grid_idx = grid_idx + 1
      ef_grid(grid_idx) = r

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
  end subroutine get_eigenfunction

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

    if (geometry == 'Cartesian') then
      eps = 1.0d0
    else if (geometry == 'cylindrical') then
      eps = r
    else
      write(*, *) "Geometry not defined correctly"
      write(*, *) "Currently set on: ", geometry
      stop
    end if

    !! \note: ic = (0.0d0, 1.0d0) and ir = (1.0d0, 0.0d0)
    select case(var)
      case('rho')
        fact = (1.0d0 / eps) * ir
      case('v1')
        fact = 1.0d0 / (ic * eps)
      case('v2')
        fact = 1.0d0 * ir
      case('v3')
        fact = (1.0d0 / eps) * ir
      case('T')
        fact = (1.0d0 / eps) * ir
      case('a1')
        fact = 1.0d0 / ic
      case('a2')
        fact = (1.0d0 / eps) * ir
      case('a3')
        fact = 1.0d0 * ir
      case default
        write(*, *) "Var not defined correctly"
        stop
    end select

  end subroutine get_transform_factor

  !> Cleaning routine, deallocates the variables in this module.
  subroutine eigenfunctions_clean()
    deallocate(ef_grid)
  end subroutine eigenfunctions_clean


end module mod_eigenfunctions
