module mod_eigenfunctions
  use mod_global_variables, only: dp, matrix_gridpts, ef_gridpts, nb_eqs, str_len_arr
  use mod_types, only: ef_type
  implicit none

  private

  real(dp), allocatable         :: ef_grid(:)
  character(str_len_arr)        :: ef_names(nb_eqs)
  type (ef_type)                :: ef_array(nb_eqs)

  public :: ef_grid, ef_names, ef_array
  public :: initialise_eigenfunctions
  public :: calculate_eigenfunctions
  public :: eigenfunctions_clean


contains

  !> Initialises the eigenfunction module, allocates X_grid.
  subroutine initialise_eigenfunctions()
    use mod_global_variables, only: write_eigenfunctions
    integer    :: i

    ! if eigenfunctions are not written to file no need to allocate them
    if (write_eigenfunctions) then
      allocate(ef_grid(ef_gridpts))
      ef_grid = 0.0d0
      do i = 1, nb_eqs
        allocate(ef_array(i) % eigenfunctions(ef_gridpts, matrix_gridpts))
      end do
    end if

    ef_names = [character(len=str_len_arr) :: 'rho', 'v1', 'v2', 'v3', 'T', 'a1', 'a2', 'a3']
  end subroutine initialise_eigenfunctions


  subroutine calculate_eigenfunctions(vr)
    complex(dp), intent(in) :: vr(matrix_gridpts, matrix_gridpts)
    complex(dp)             :: ef_values(ef_gridpts)
    integer                 :: i, j

    do j = 1, nb_eqs
      ! eigenfunction index
      ef_array(j) % index = j
      ! eigenfunction name corresponding to the index
      ef_array(j) % name = ef_names(j)

      do i = 1, matrix_gridpts
        ! the eigenfunction for each eigenvalue is stored in the column of 'eigenfunctions'
        ! with the same index as omega, so column indices correspond to the eigenvalue at that index.
        call get_eigenfunction(trim(ef_array(j) % name), vr(:, i), ef_values)
        ef_array(j) % eigenfunctions(:, i) = ef_values
      end do
    end do
  end subroutine calculate_eigenfunctions

  !> Subroutine to calculate the eigenfunction of a given variable, taking
  !! the finite elements into account.
  !! @param[in] var   The variable (rho, v1, ...) for the eigenfunction
  !! @param[in] eigenvector The eigenvector for which to calculate the
  !!                        eigenfunction. This is one of the columns in 'vr'
  !! @param[out] Y  The eigenfunction corresponding to the requested variable
  !!                and eigenvector
  subroutine get_eigenfunction(name, eigenvector, Y)
    use mod_global_variables, only: dim_subblock, gridpts
    use mod_grid, only: grid

    character(len=*), intent(in)  :: name
    complex(dp), intent(in)       :: eigenvector(matrix_gridpts)
    complex(dp), intent(out)      :: Y(ef_gridpts)

    integer                       :: idx1, idx2, grid_idx, i
    real(dp)                      :: r, r_lo, r_hi
    real(dp)                      :: spline(4)
    complex(dp)                   :: fact

    Y = (0.0d0, 0.0d0)

    select case(name)
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
        write(*, *) "Currently set on: ", name
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

    call get_transform_factor(name, r, fact)
    call get_spline(name, r, r_lo, r_hi, spline)
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

      call get_transform_factor(name, r, fact)
      call get_spline(name, r, r_lo, r_hi, spline)
      Y(grid_idx) = Y(grid_idx) + fact*eigenvector(idx1)*spline(2) &
                                + fact*eigenvector(idx2)*spline(4)
      Y(grid_idx) = Y(grid_idx) + fact*eigenvector(idx1+dim_subblock)*spline(1) &
                                + fact*eigenvector(idx2+dim_subblock)*spline(3)

      ! Contribution from end point
      r = r_hi
      grid_idx = grid_idx + 1
      ef_grid(grid_idx) = r

      call get_transform_factor(name, r, fact)
      call get_spline(name, r, r_lo, r_hi, spline)
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
  subroutine get_spline(name, r, r_lo, r_hi, spline)
    use mod_spline_functions, only: cubic_factors, quadratic_factors

    character(len=*), intent(in)    :: name
    real(dp), intent(in)            :: r, r_lo, r_hi
    real(dp), intent(out)           :: spline(4)

    if ((name == "v1") .or. (name == "a2") .or. (name == "a3")) then
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
  subroutine get_transform_factor(name, r, fact)
    use mod_global_variables, only: geometry, ir, ic

    character(len=*), intent(in) :: name
    real(dp), intent(in)         :: r
    complex(dp), intent(out)     :: fact
    real(dp)                     :: eps

    if (geometry == 'Cartesian') then
      eps = 1.0d0
    else
      eps = r
    end if

    !! \note: ic = (0.0d0, 1.0d0) and ir = (1.0d0, 0.0d0)
    select case(name)
      case('rho')
        fact = ir / eps
      case('v1')
        fact = - ic / eps
      case('v2')
        fact = ir
      case('v3')
        fact = ir / eps
      case('T')
        fact = ir / eps
      case('a1')
        fact = - ic
      case('a2')
        fact = ir / eps
      case('a3')
        fact = ir
      case default
        write(*, *) "Eigenfunction name not defined correctly"
        stop
    end select

  end subroutine get_transform_factor

  !> Cleaning routine, deallocates the variables in this module.
  subroutine eigenfunctions_clean()
    use mod_global_variables, only: write_eigenfunctions

    integer   :: i

    if (write_eigenfunctions) then
      deallocate(ef_grid)
      do i = 1, nb_eqs
        deallocate(ef_array(i) % eigenfunctions)
      end do
    end if
  end subroutine eigenfunctions_clean


end module mod_eigenfunctions
