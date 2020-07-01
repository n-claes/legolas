! =============================================================================
!> @brief   Module containing everything eigenfunction-related.
!! @details This module composes the eigenfunctions based on the eigenvectors
!!          coming out of the solver routines. Assembling the eigenfunctions is needed
!!          because of the finite element matrix representation.
!!          In every interval of the base grid (i, i+1) eigenfunctions are saved at the
!!          left and right edges, as well as in the centre of the interval.
!!          This results in eigenfunction arrays with a size <em> 2*gridpts - 1 </em>.
module mod_eigenfunctions
  use mod_global_variables, only: dp, matrix_gridpts, ef_gridpts, nb_eqs, str_len_arr
  use mod_types, only: ef_type
  implicit none

  private

  !> grid to which the eigenfunction is assembled
  real(dp), allocatable         :: ef_grid(:)
  !> array containing the eigenfunction names as strings
  character(str_len_arr)        :: ef_names(nb_eqs)
  !> type containig all eigenfunctions
  type (ef_type)                :: ef_array(nb_eqs)

  public :: ef_grid, ef_names, ef_array
  public :: initialise_eigenfunctions
  public :: calculate_eigenfunctions
  public :: eigenfunctions_clean

contains


  !> @brief   Main initialisations of this module.
  !! @details Allocates and initialises the eigenfunction grid, types and names.
  !! @note    This is only done if the global variable for eigenfunction assembly is \p True.
  !!          If eigenfunctions are not saved anyway it makes no sense to assemble them.
  subroutine initialise_eigenfunctions()
    use mod_global_variables, only: write_eigenfunctions
    integer    :: i

    if (write_eigenfunctions) then
      allocate(ef_grid(ef_gridpts))
      ef_grid = 0.0d0
      do i = 1, nb_eqs
        allocate(ef_array(i) % eigenfunctions(ef_gridpts, matrix_gridpts))
      end do
    end if

    ef_names = [character(len=str_len_arr) :: 'rho', 'v1', 'v2', 'v3', 'T', 'a1', 'a2', 'a3']
  end subroutine initialise_eigenfunctions


  !> @brief   Calculates the eigenfunctions.
  !! @details Assembles the eigenfunctions for every eigenvalue and variable,
  !!          based on the right eigenvectors. The eigenfunctions are transformed
  !!          back to their 'actual' values
  !! @note    All eigenfunctions with a value smaller than <tt>1e-10</tt> are set to zero.
  !!          This is separately checked for the real and imaginary parts.
  !! @param[in] vr  the matrix of right eigenvectors
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


  !> @brief   Calculates a single eigenfunction based on an eigenvector.
  !! @details Assembles the eigenfunction corresponding to \p ef_idx and the
  !!          given eigenvector. The eigenvector supplied is the one
  !!          corresponding to one particular eigenvalue.
  !! @param[in] ef_idx        the index of the variable in the state vector
  !! @param[in] eigenvector   the eigenvector for this particular eigenvalue
  !! @param[out] Y            the assembled eigenfunction
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


  !> @brief   Retrieves the spline for a particular variable.
  !! @details Calls the correct finite element basis function depending on the
  !!          variable passed.
  !! @param[in] ef_idx  the index of the variable in the state vector
  !! @param[in] r       the current position in the grid
  !! @param[in] r_lo    left edge of the current grid interval
  !! @param[in] r_hi    right edge of the current grid interval
  !! @param[out] spline the finite element basis function, depending on \p ef_idx
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


  !! @brief   Re-transformation of the eigenfunction.
  !! @details Transforms the eigenfunction back to its 'actual' value.
  !!          For example, \p rho is defined as <tt>eps*rho</tt> in the equations,
  !!          in this case the eigenfunction is divided by the scale factor.
  !! @exception  Error  If an unknown ef_idx is supplied.
  !! @param[in] ef_idx  index of the variable in the state vector
  !! @param[in, out] ef_values   eigenfunction values, transformed on return
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


  !> @brief   Module cleaning routine.
  !! @details Deallocates the eigenfunction arrays.
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
