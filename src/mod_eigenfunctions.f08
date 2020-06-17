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
    use mod_check_values, only: check_small_values

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
        call get_eigenfunction(ef_array(j) % index, vr(:, i), ef_values)
        ! undo variable transformation for 'actual' eigenfunction
        call transform_eigenfunction(ef_array(j) % index, ef_values)
        call check_small_values(ef_values, tol=1.0d-10)
        ef_array(j) % eigenfunctions(:, i) = ef_values
      end do
    end do
  end subroutine calculate_eigenfunctions

  !> Subroutine to calculate the eigenfunction of a given variable, taking
  !! the finite elements into account.
  !! @param[in] ef_idx   The index of the eigenfunction in the ef_names array
  !! @param[in] eigenvector The eigenvector for which to calculate the
  !!                        eigenfunction. This is one of the columns in 'vr'
  !! @param[out] Y  The eigenfunction corresponding to the requested variable
  !!                and eigenvector
  subroutine get_eigenfunction(ef_idx, eigenvector, Y)
    use mod_global_variables, only: dim_subblock, gridpts
    use mod_grid, only: grid
    use mod_logging, only: log_message

    integer, intent(in)           :: ef_idx
    complex(dp), intent(in)       :: eigenvector(matrix_gridpts)
    complex(dp), intent(out)      :: Y(ef_gridpts)

    integer                       :: idx1, idx2, grid_idx, i
    real(dp)                      :: r, r_lo, r_hi
    real(dp)                      :: spline(4)

    ! initialise eigenfunction to zero
    Y = (0.0d0, 0.0d0)
    ! map ef_idx to actual subblock index. So 1 -> 1, 2 -> 3, 3 -> 5 etc.
    idx1 = 2 * ef_idx - 1
    idx2 = idx1 + 1

    !! note: eigenfunctions are saved in each interval [grid(i), grid(i+1)] at
    !! the left edge, centre and right edge

    ! Contribution from first gridpoint, left edge
    r_lo = grid(1)
    r_hi = grid(2)
    r    = grid(1)
    grid_idx = 1
    ef_grid(grid_idx) = r
    call get_spline(ef_idx, r, r_lo, r_hi, spline)
    Y(grid_idx) = Y(grid_idx) + eigenvector(idx1) * spline(2) &
                              + eigenvector(idx2) * spline(4) &
                              + eigenvector(idx1 + dim_subblock) * spline(1) &
                              + eigenvector(idx2 + dim_subblock) * spline(3)

    ! Contribution from other gridpoints
    do i = 1, gridpts-1
      ! Contribution from centre
      r_lo = grid(i)
      r_hi = grid(i + 1)
      r    = 0.5 * (r_lo + r_hi)
      ! save gridpoint
      grid_idx = grid_idx + 1
      ef_grid(grid_idx) = r

      call get_spline(ef_idx, r, r_lo, r_hi, spline)
      Y(grid_idx) = Y(grid_idx) + eigenvector(idx1) * spline(2) &
                                + eigenvector(idx2) * spline(4) &
                                + eigenvector(idx1 + dim_subblock) * spline(1) &
                                + eigenvector(idx2 + dim_subblock) * spline(3)

      ! Contribution from end point
      r = r_hi
      grid_idx = grid_idx + 1
      ef_grid(grid_idx) = r
      call get_spline(ef_idx, r, r_lo, r_hi, spline)
      Y(grid_idx) = Y(grid_idx) + eigenvector(idx1) * spline(2) &
                                + eigenvector(idx2) * spline(4) &
                                + eigenvector(idx1 + dim_subblock) * spline(1) &
                                + eigenvector(idx2 + dim_subblock) * spline(3)

      ! Increment indices of eigenvector elements
      idx1 = idx1 + dim_subblock
      idx2 = idx2 + dim_subblock
    end do
  end subroutine get_eigenfunction

  !> Calls the right finite element spline depending on the variable passed.
  !! @param[in] ef_idx   The index of the eigenfunction in the ef_names array
  !! @param[in] r     The current position in the grid
  !! @param[in] r_lo  Left edge of the current grid interval
  !! @param[in] r_hi  Right edge of the current grid interval
  !! @param[out] spline The finite element basis function, depending on 'var'
  subroutine get_spline(ef_idx, r, r_lo, r_hi, spline)
    use mod_spline_functions, only: cubic_factors, quadratic_factors

    integer, intent(in)             :: ef_idx
    real(dp), intent(in)            :: r, r_lo, r_hi
    real(dp), intent(out)           :: spline(4)

    if (ef_idx == 2 .or. ef_idx == 7 .or. ef_idx == 8) then
      call cubic_factors(r, r_lo, r_hi, spline)
    else
      call quadratic_factors(r, r_lo, r_hi, spline)
    end if
  end subroutine get_spline


  !> Subroutine to retrieve the transformation factor. For example, rho was
  !! rescaled by eps*rho = RHO, such that in this case 'fact' will be 1/eps.
  !! @param[in] ef_idx   index of the eigenfunction name in the ef_names array
  !! @param[in, out] ef_values    the eigenfunction values, transformed on exit
  subroutine transform_eigenfunction(ef_idx, ef_values)
    use mod_global_variables, only: ic, geometry
    use mod_logging, only: log_message

    integer, intent(in)         :: ef_idx
    complex(dp), intent(inout)  :: ef_values(ef_gridpts)
    real(dp)     :: ef_eps(ef_gridpts)

    ! we can't use eps_grid here from mod_grid, since the eigenfunctions
    ! were interpolated between gridpoints as well
    if (geometry == 'Cartesian') then
      ef_eps = 1.0d0
    else
      ef_eps = ef_grid
    end if

    select case(ef_idx)
    case(1) ! rho -> eps*rho
      ef_values = ef_values / ef_eps
    case(2) ! v1 -> i*eps*v1
      ef_values = ef_values / (ef_eps * ic)
    case(3) ! v2 -> v2
      ! do nothing
    case(4) ! v3 -> eps*v3
      ef_values = ef_values / ef_eps
    case(5) ! T1 -> eps*T1
      ef_values = ef_values / ef_eps
    case(6) ! a1 -> i*a1
      ef_values = ef_values / ic
    case(7) ! a2 -> eps*a2
      ef_values = ef_values / ef_eps
    case(8) ! a3 -> a3
      ! do nothing
    case default
      call log_message('wrong eigenfunction index in transform_eigenfunction()', level='error')
    end select
  end subroutine transform_eigenfunction

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
