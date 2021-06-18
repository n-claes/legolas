! =============================================================================
!> This module composes perturbed quantities of interest based on the
!! eigenfunctions and eigenvectors. As a result, all quantities are defined on
!! the same grid as the eigenfunctions, ef_grid.
!! The computed quantities are
!!    - div v1
!!    - all 3 components of curl v1
!!    - all 3 components of B1
!!    - div B1
!!    - all 3 components of curl B1

module mod_postprocessing
  use mod_global_variables, only: dp, ic, matrix_gridpts, ef_gridpts, str_len_arr
  use mod_equilibrium_params, only: k2, k3
  use mod_eigenfunctions, only: ef_grid, ef_array
  use mod_types, only: pp_type
  implicit none

  private

  !> number of postprocessing quantities
  integer                       :: nb_pp = 11
  !> array containing the postprocessed quantity names as strings
  character(str_len_arr)        :: pp_names(11)
  !> type containig all postprocessed quantities
  type (pp_type)                :: pp_array(11)

  public :: pp_names, pp_array
  public :: initialise_postprocessing
  public :: calculate_postprocessed
  public :: postprocessing_clean

contains


  !> Main initialisations of this module.
  !! Allocates and initialises the types and names.
  !! Passing the optional argument <tt>nb_evs</tt> sets the number of eigenvalues
  !! that are calculated and limits the size of the eigenfunction arrays accordingly.
  !! This routine is only called if <tt>write_postprocessed = .true.</tt>.
  subroutine initialise_postprocessing(nb_evs)
    !> the number of eigenvalues that are calculated, defaults to all (matrix dim)
    integer, intent(in), optional :: nb_evs

    integer    :: i, nev

    if (present(nb_evs)) then
      nev = nb_evs
    else
      nev = matrix_gridpts
    end if

    do i = 1, nb_pp
      allocate(pp_array(i) % quantities(ef_gridpts, nev))
    end do

    pp_names = [ &
      character(len=str_len_arr) :: 'div v1', '(curl v1)_1', '(curl v1)_2', &
                                    '(curl v1)_3', 'B11', 'B12', 'B13', 'div B1', &
                                    '(curl B1)_1', '(curl B1)_2', '(curl B1)_3' &
    ]
  end subroutine initialise_postprocessing


  !> Calculates the postprocessed quantities for every eigenvalue.
  subroutine calculate_postprocessed(vr)
    use mod_check_values, only: check_small_values

    !> the matrix of right eigenvectors
    complex(dp), intent(in) :: vr(:, :)
    !> quantities values
    complex(dp) :: single(ef_gridpts), triple(ef_gridpts, 3)
    integer     :: i, j

    do j = 1, nb_pp
      ! quantity index
      pp_array(j) % index = j
      ! quantity name corresponding to the index
      pp_array(j) % name = pp_names(j)
    end do

    !! Divergence of v1
    do i = 1, size(vr, dim=2)
      call get_div_velocity(i, vr(:, i), single)
      pp_array(1) % quantities(:, i) = single
    end do

    !! Curl of v1
    do i = 1, size(vr, dim=2)
      call get_curl_velocity(i, vr(:, i), triple)
      pp_array(2) % quantities(:, i) = triple(:, 1)
      pp_array(3) % quantities(:, i) = triple(:, 2)
      pp_array(4) % quantities(:, i) = triple(:, 3)
    end do

    !! Magnetic field B1
    do i = 1, size(vr, dim=2)
      call get_magnetic_field(i, vr(:, i), triple)
      pp_array(5) % quantities(:, i) = triple(:, 1)
      pp_array(6) % quantities(:, i) = triple(:, 2)
      pp_array(7) % quantities(:, i) = triple(:, 3)
    end do

    !! Divergence of B1
    do i = 1, size(vr, dim=2)
      call get_div_magnetic(i, vr(:, i), single)
      pp_array(8) % quantities(:, i) = single
    end do

    !! Curl of B1
    do i = 1, size(vr, dim=2)
      call get_curl_magnetic(i, vr(:, i), triple)
      pp_array(9) % quantities(:, i) = triple(:, 1)
      pp_array(10) % quantities(:, i) = triple(:, 2)
      pp_array(11) % quantities(:, i) = triple(:, 3)
    end do
  end subroutine calculate_postprocessed


  !> Calculates the divergence of the perturbed velocity v1
  subroutine get_div_velocity(index, eigenvector, divergence)
    integer, intent(in)       :: index
    complex(dp), intent(in)   :: eigenvector(matrix_gridpts)
    complex(dp), intent(out)  :: divergence(ef_gridpts)

    !! Derivatives are transformed variables, regular variables are 'actual' values
    complex(dp) :: v2(ef_gridpts), v3(ef_gridpts)
    complex(dp) :: dv1(ef_gridpts)
    integer     :: m
    real(dp)    :: eps, deps

    v2 = ef_array(3) % eigenfunctions(:, index)
    v3 = ef_array(4) % eigenfunctions(:, index)
    call get_eigenfunction_deriv(2, eigenvector, dv1, 1)

    do m = 1, ef_gridpts
      call set_eps(m, eps, deps)
      divergence(m) = ic * (-dv1(m) / eps + k2 * v2(m) / eps + k3 * v3(m))
    end do
  end subroutine get_div_velocity


  !> Calculates the curl of the perturbed velocity v1
  subroutine get_curl_velocity(index, eigenvector, curl)
    integer, intent(in)       :: index
    complex(dp), intent(in)   :: eigenvector(matrix_gridpts)
    complex(dp), intent(out)  :: curl(ef_gridpts, 3)

    !! Derivatives are transformed variables, regular variables are 'actual' values
    complex(dp) :: v1(ef_gridpts), v2(ef_gridpts), v3(ef_gridpts)
    complex(dp) :: dv2(ef_gridpts), dv3(ef_gridpts)
    integer     :: m
    real(dp)    :: eps, deps

    v1 = ef_array(2) % eigenfunctions(:, index)
    v2 = ef_array(3) % eigenfunctions(:, index)
    v3 = ef_array(4) % eigenfunctions(:, index)
    call get_eigenfunction_deriv(3, eigenvector, dv2, 1)
    call get_eigenfunction_deriv(4, eigenvector, dv3, 1)

    do m = 1, ef_gridpts
      call set_eps(m, eps, deps)
      curl(m, 1) = ic * (k2 * v3(m) / eps - k3 * v2(m))
      curl(m, 2) = ic * k3 * v1(m) - dv3(m) / eps - deps * v3(m) / eps
      curl(m, 3) = dv2(m) + (deps * v2(m) - ic * k2 * v1(m)) / eps
    end do
  end subroutine get_curl_velocity


  !> Calculates the perturbed magnetic field B1
  subroutine get_magnetic_field(index, eigenvector, field)
    integer, intent(in)       :: index
    complex(dp), intent(in)   :: eigenvector(matrix_gridpts)
    complex(dp), intent(out)  :: field(ef_gridpts, 3)

    !! Derivatives are transformed variables, regular variables are 'actual' values
    complex(dp) :: a1(ef_gridpts), a2(ef_gridpts), a3(ef_gridpts)
    complex(dp) :: da2(ef_gridpts), da3(ef_gridpts)
    integer     :: m
    real(dp)    :: eps, deps

    a1 = ef_array(6) % eigenfunctions(:, index)
    a2 = ef_array(7) % eigenfunctions(:, index)
    a3 = ef_array(8) % eigenfunctions(:, index)
    call get_eigenfunction_deriv(7, eigenvector, da2, 1)
    call get_eigenfunction_deriv(8, eigenvector, da3, 1)

    do m = 1, ef_gridpts
      call set_eps(m, eps, deps)
      field(m, 1) = ic * (k2 * a3(m) / eps - k3 * a2(m))
      field(m, 2) = ic * k3 * a1(m) - da3(m)
      field(m, 3) = (da2(m) - ic * k2 * a1(m)) / eps
    end do
  end subroutine get_magnetic_field


  !> Calculates the divergence of the perturbed magnetic field B1
  subroutine get_div_magnetic(index, eigenvector, divergence)
    integer, intent(in)       :: index
    complex(dp), intent(in)   :: eigenvector(matrix_gridpts)
    complex(dp), intent(out)  :: divergence(ef_gridpts)

    !! Derivatives are transformed variables, regular variables are 'actual' values
    complex(dp) :: a1(ef_gridpts), a2(ef_gridpts), a3(ef_gridpts)
    complex(dp) :: da2(ef_gridpts), da3(ef_gridpts)
    integer     :: m
    real(dp)    :: eps, deps

    a1 = ef_array(6) % eigenfunctions(:, index)
    a2 = ef_array(7) % eigenfunctions(:, index)
    a3 = ef_array(8) % eigenfunctions(:, index)
    call get_eigenfunction_deriv(7, eigenvector, da2, 1)
    call get_eigenfunction_deriv(8, eigenvector, da3, 1)

    do m = 1, ef_gridpts
      call set_eps(m, eps, deps)
      divergence(m) = ic * (k2 * da3(m) - k3 * da2(m)) / eps &
                      + ic * k2 * (ic * k3 * a1(m) - da3(m)) / eps &
                      + ic * k3 * (da2(m) - ic * k2 * a1(m)) / eps
    end do
  end subroutine get_div_magnetic


  !> Calculates the curl of the perturbed magnetic field B1
  subroutine get_curl_magnetic(index, eigenvector, curl)
    integer, intent(in)       :: index
    complex(dp), intent(in)   :: eigenvector(matrix_gridpts)
    complex(dp), intent(out)  :: curl(ef_gridpts, 3)

    !! Derivatives are transformed variables, regular variables are 'actual' values
    complex(dp) :: a1(ef_gridpts), a2(ef_gridpts), a3(ef_gridpts)
    complex(dp) :: da1(ef_gridpts), da2(ef_gridpts), da3(ef_gridpts)
    complex(dp) :: dda2(ef_gridpts), dda3(ef_gridpts)
    integer     :: m
    real(dp)    :: eps, deps

    a1 = ef_array(6) % eigenfunctions(:, index)
    a2 = ef_array(7) % eigenfunctions(:, index)
    a3 = ef_array(8) % eigenfunctions(:, index)
    call get_eigenfunction_deriv(6, eigenvector, da1, 1)
    call get_eigenfunction_deriv(7, eigenvector, da2, 1)
    call get_eigenfunction_deriv(8, eigenvector, da3, 1)
    call get_eigenfunction_deriv(7, eigenvector, dda2, 2)
    call get_eigenfunction_deriv(8, eigenvector, dda3, 2)

    do m = 1, ef_gridpts
      call set_eps(m, eps, deps)
      curl(m, 1) = ic * k2 * (da2(m) - ic * k2 * a1(m)) / eps**2 &
                   - ic * k3 * (ic * k3 * a1(m) - da3(m))
      curl(m, 2) = -k3 * (k2 * a3(m) / eps - k3 * a2(m)) &
                   - (dda2(m) - k2 * da1(m)) / eps &
                   + deps * (da2(m) - ic * k2 * a1(m)) / eps**2
      curl(m, 3) = k3 * da1(m) - dda3(m) + deps * (ic * k3 * a1(m) - da3(m)) / eps &
                   + k2 * (k2 * a3(m) / eps - k3 * a2(m)) / eps
    end do
  end subroutine get_curl_magnetic


  !> Calculates the derivative of a single eigenfunction based on an eigenvector.
  !! The eigenvector supplied is the one corresponding to one particular
  !! eigenvalue.
  subroutine get_eigenfunction_deriv(ef_idx, eigenvector, Y, deriv)
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
    !> the deriv-th spline derivative, defaults to 0
    integer, optional             :: deriv

    if (.not. present(deriv)) then
      deriv = 0
    end if

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
    call get_spline_deriv(ef_idx, r, r_lo, r_hi, spline, deriv)
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

      call get_spline_deriv(ef_idx, r, r_lo, r_hi, spline, deriv)
      Y(grid_idx) = Y(grid_idx) + eigenvector(idx1) * spline(2) &
                                + eigenvector(idx2) * spline(4) &
                                + eigenvector(idx1 + dim_subblock) * spline(1) &
                                + eigenvector(idx2 + dim_subblock) * spline(3)

      ! Contribution from end point
      r = r_hi
      grid_idx = grid_idx + 1
      ef_grid(grid_idx) = r
      call get_spline_deriv(ef_idx, r, r_lo, r_hi, spline, deriv)
      Y(grid_idx) = Y(grid_idx) + eigenvector(idx1) * spline(2) &
                                + eigenvector(idx2) * spline(4) &
                                + eigenvector(idx1 + dim_subblock) * spline(1) &
                                + eigenvector(idx2 + dim_subblock) * spline(3)

      ! Increment indices of eigenvector elements
      idx1 = idx1 + dim_subblock
      idx2 = idx2 + dim_subblock
    end do
  end subroutine get_eigenfunction_deriv


  !> Retrieves the spline or one of its derivatives for a particular variable by
  !! calling the correct finite element basis function depending on the
  !! variable passed.
  subroutine get_spline_deriv(ef_idx, r, r_lo, r_hi, spline, deriv)
    use mod_spline_functions, only: cubic_factors, quadratic_factors, &
                                    cubic_factors_deriv, quadratic_factors_deriv, &
                                    cubic_factors_deriv2

    !> the index of the variable in the state vector
    integer, intent(in)     :: ef_idx
    !> the current position in the grid
    real(dp), intent(in)    :: r
    !> left edge of the current grid interval
    real(dp), intent(in)    :: r_lo
    !> right edge of the current grid interval
    real(dp), intent(in)    :: r_hi
    !> the finite element basis functions for this grid position
    real(dp), intent(out)   :: spline(4)
    !> the deriv-th spline derivative, defaults to 0 (no derivative)
    integer, optional       :: deriv

    if (.not. present(deriv)) then
      deriv = 0
    end if

    if (ef_idx == 2 .or. ef_idx == 7 .or. ef_idx == 8) then
      if (deriv == 2) then
        call cubic_factors_deriv2(r, r_lo, r_hi, spline)
      else if (deriv == 1) then
        call cubic_factors_deriv(r, r_lo, r_hi, spline)
      else
        call cubic_factors(r, r_lo, r_hi, spline)
      end if
    else
      if (deriv == 1) then
        call quadratic_factors_deriv(r, r_lo, r_hi, spline)
      else
        call quadratic_factors(r, r_lo, r_hi, spline)
      end if
    end if
  end subroutine get_spline_deriv


  subroutine set_eps(index, eps, deps)
    use mod_global_variables, only: geometry

    integer, intent(in)   :: index
    real(dp), intent(out) :: eps, deps

    if (geometry == 'cylindrical') then
      eps  = ef_grid(index)
      deps = 1.0d0
    else
      eps  = 1.0d0
      deps = 1.0d0
    end if
  end subroutine set_eps


  !> Cleaning routine, deallocates the postprocessing arrays.
  subroutine postprocessing_clean()
    integer   :: i

    if (allocated(pp_array(1) % quantities)) then
      do i = 1, nb_pp
        deallocate(pp_array(i) % quantities)
      end do
    end if
  end subroutine postprocessing_clean

end module mod_postprocessing
