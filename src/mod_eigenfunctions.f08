! =============================================================================
!> This module composes the eigenfunctions based on the eigenvectors
!! coming out of the solver routines. Assembling the eigenfunctions is needed
!! because of the finite element matrix representation.
!! In every interval of the base grid \((i, i+1)\) eigenfunctions are saved at the
!! left and right edges, as well as in the centre of the interval.
!! This results in eigenfunction arrays with a size \(2*gridpts - 1\).
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


  !> Main initialisations of this module.
  !! Allocates and initialises the eigenfunction grid, types and names.
  !! Passing the optional argument <tt>nb_evs</tt> sets the number of eigenvalues
  !! that are calculated and limits the size of the eigenfunction arrays accordingly.
  !! This routine is only called if <tt>write_eigenfunctions = .true.</tt>.
  subroutine initialise_eigenfunctions(nb_evs)
    !> the number of eigenvalues that are calculated, defaults to all (matrix dim)
    integer, intent(in), optional :: nb_evs

    integer    :: i, nev

    if (present(nb_evs)) then
      nev = nb_evs
    else
      nev = matrix_gridpts
    end if

    allocate(ef_grid(ef_gridpts))
    ef_grid = 0.0d0
    do i = 1, nb_eqs
      allocate(ef_array(i) % eigenfunctions(ef_gridpts, nev))
    end do

    ef_names = [ &
      character(len=str_len_arr) :: 'rho', 'v1', 'v2', 'v3', 'T', 'a1', 'a2', 'a3' &
    ]
  end subroutine initialise_eigenfunctions


  !> Calculates the eigenfunctions for every eigenvalue and variable,
  !! based on the right eigenvectors. The eigenfunctions are transformed
  !! back to their 'actual' values
  subroutine calculate_eigenfunctions(vr)
    !> the matrix of right eigenvectors
    complex(dp), intent(in) :: vr(:, :)
    !> eigenfunctio values
    complex(dp)             :: ef_values(ef_gridpts)
    integer                 :: i, j

    do j = 1, nb_eqs
      ! eigenfunction index
      ef_array(j) % index = j
      ! eigenfunction name corresponding to the index
      ef_array(j) % name = ef_names(j)

      do i = 1, size(vr, dim=2)
        ! the eigenfunction for each eigenvalue is stored in the column of
        ! 'eigenfunctions' with the same index as omega,
        ! so column indices correspond to the eigenvalue at that index.
        call get_eigenfunction(ef_array(j) % index, vr(:, i), ef_values)
        ! undo variable transformation for 'actual' eigenfunction
        call transform_eigenfunction(ef_array(j) % index, ef_values)
        ef_array(j) % eigenfunctions(:, i) = ef_values
      end do
    end do
  end subroutine calculate_eigenfunctions


  !> Calculates a single eigenfunction based on an eigenvector.
  !! Assembles the eigenfunction corresponding to <tt>ef_idx</tt> and the
  !! given eigenvector. The eigenvector supplied is the one
  !! corresponding to one particular eigenvalue.
  subroutine get_eigenfunction(ef_idx, eigenvector, Y)
    use mod_global_variables, only: dim_subblock, gridpts
    use mod_grid, only: grid
    use mod_logging, only: log_message

    !> the index of the variable in the state vector
    integer, intent(in)           :: ef_idx
    !> the eigenvector for this particular eigenvalue
    complex(dp), intent(in)       :: eigenvector(matrix_gridpts)
    !> the assembled eigenfunction
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


  !> Retrieves the spline for a particular variable by
  !! calling the correct finite element basis function depending on the
  !! variable passed.
  subroutine get_spline(ef_idx, r, r_lo, r_hi, spline)
    use mod_spline_functions, only: cubic_factors, quadratic_factors

    !> the index of the variable in the state vector
    integer, intent(in)             :: ef_idx
    !> the current position in the grid
    real(dp), intent(in)            :: r
    !> left edge of the current grid interval
    real(dp), intent(in)            :: r_lo
    !> right edge of the current grid interval
    real(dp), intent(in)            :: r_hi
    !> the finite element basis functions for this grid position
    real(dp), intent(out)           :: spline(4)

    if (ef_idx == 2 .or. ef_idx == 7 .or. ef_idx == 8) then
      call cubic_factors(r, r_lo, r_hi, spline)
    else
      call quadratic_factors(r, r_lo, r_hi, spline)
    end if
  end subroutine get_spline


  !> Re-transformation of the eigenfunction.
  !! Transforms the eigenfunction back to its 'actual' value.
  !! For example, \(\rho\) is defined as \(\varepsilon\rho\) in the equations,
  !! in this case the eigenfunction is divided by the scale factor.
  !! @warning Throws an error if an unknown <tt>ef_idx</tt> is supplied.
  subroutine transform_eigenfunction(ef_idx, ef_values)
    use mod_global_variables, only: ic, geometry
    use mod_logging, only: log_message

    !> index of the variable in the state vector
    integer, intent(in)         :: ef_idx
    !> eigenfunction values, transformed on exit
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


  !> Cleaning routine, deallocates the eigenfunction arrays.
  subroutine eigenfunctions_clean()
    integer   :: i

    if (allocated(ef_grid)) then
      deallocate(ef_grid)
      do i = 1, nb_eqs
        deallocate(ef_array(i) % eigenfunctions)
      end do
    end if
  end subroutine eigenfunctions_clean

end module mod_eigenfunctions
