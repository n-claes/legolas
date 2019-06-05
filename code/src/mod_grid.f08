module mod_grid
  use mod_global_variables
  implicit none

  real(dp), allocatable      :: grid(:), grid_gauss(:)

contains

  subroutine initialise_grid()
    integer                  :: i
    real(dp)                 :: dx

    allocate(grid(gridpts))
    allocate(grid_gauss(4*gridpts))

    ! Initialise grids
    grid       = 0.0d0
    grid_gauss = 0.0d0

    ! minus one here to include x_end
    dx = (x_end - x_start) / (gridpts-1)
    do i = 1, gridpts
      grid(i) = x_start + (i - 1)*dx
    end do

    if (mesh_accumulation) then
      call accumulate_mesh()
    end if

  end subroutine initialise_grid

  !> Subroutine to re-grid the mesh to a non-uniform spacing.
  !! This is based on two Gaussian curves with known widths (from
  !! mod_global_variables); the grid is accumulated near each Gaussian maximum,
  !! using equidistribution based on the integral under the curve defined
  !! by the function gaussian().
  subroutine accumulate_mesh()
    integer                  :: i, integral_gridpts_1, integral_gridpts_2
    real(dp)                 :: dx, dx_0, xi, bgf, fact, dx_eq
    real(dp)                 :: gauss_xi, gauss_xi_eq
    real(dp)                 :: x_sum, x_sum_prev, x_norm
    real(dp)                 :: xi_weighted

    print*,"Redefining grid with mesh accumulation"

    if (integral_gridpts /= gridpts - 1) then
      stop "WARNING: integral gridpoints must be gridpoints - 1."
    end if

    bgf  = 0.3d0 !background
    fact = 1.0d0

    ! first evaluation of integral to get weighted values
    integral_gridpts_1 = 2*integral_gridpts + 1
    dx = (grid(gridpts) - grid(1)) / float(integral_gridpts_1 - 1)
    xi = grid(1)
    x_sum = 0.0d0

    do i = 1, integral_gridpts_1
      gauss_xi = gaussian(xi, bgf, fact)
      x_sum   = x_sum + (gauss_xi * dx)
      xi      = xi + dx
    end do
    x_norm = (grid(gridpts) - grid(1)) / x_sum

    ! second evaluation of integral using weighted points
    integral_gridpts_2 = 50*integral_gridpts
    dx_eq    = (grid(gridpts) - grid(1)) / float(integral_gridpts)
    xi       = grid(1)
    x_sum    = 0.0d0           ! x0 here
    gauss_xi = gaussian(xi, bgf, fact) * x_norm   ! This is at x0 for now
    dx_0     = (grid(gridpts) - grid(1)) * gauss_xi / float(integral_gridpts_2)


    do i = 2, integral_gridpts
      gauss_xi_eq = float(i - 1) * dx_eq + grid(1)
      do while (gauss_xi_eq > x_sum)  !x_sum is 0 at first pass
       dx         = dx_0 / gauss_xi
       xi         = xi + dx
       x_sum_prev = x_sum
       gauss_xi   = gaussian(xi, bgf, fact)
       gauss_xi   = gauss_xi * x_norm
       x_sum      = x_sum + (gauss_xi * dx)
      end do

      xi_weighted = (gauss_xi_eq - x_sum_prev) / (x_sum - x_sum_prev)

      ! Re-define grid
      grid(i) = xi - dx*(1.0d0 - xi_weighted)
    end do

    ! Ensure correct end points and final spacing
    grid(1) = x_start
    grid(integral_gridpts+1) = x_end
    grid(integral_gridpts) = 0.5 * (grid(integral_gridpts - 1) &
                                    + grid(integral_gridpts + 1))

  end subroutine accumulate_mesh



  !> Function to calculate a Gaussian curve based on known widths and
  !! expected values (from mod_global_variables). The Gaussian is evaluated
  !! in x.
  function gaussian(x, bgf, fact) result(f_gauss)

    real(dp), intent(in)    :: x, bgf, fact
    real(dp)                :: f_gauss
    real(dp)                :: gauss_1, gauss_2, norm_1, norm_2

    norm_1 = 1.0d0 / (sigma_1 * sqrt(2*dpi))
    norm_2 = 1.0d0 / (sigma_2 * sqrt(2*dpi))

    gauss_1 = norm_1 * exp(-0.5d0 * ((x - ev_1) / sigma_1)**2)
    gauss_2 = norm_2 * exp(-0.5d0 * ((x - ev_2) / sigma_2)**2)

    f_gauss = bgf + (1.0d0 - bgf) * (gauss_1 + fact*gauss_2) / fact

  end function gaussian


  !> Subroutine to calculate the derivative of a given array (y_array) with
  !! respect to x_array using central difference discretization.
  !! @param x_array       Used to calculate dx for dy/dx.
  !! @param y_array       Used to calculate dy for dy/dx.
  !! @param d_y_array_dx  y_array derived with respect to x_array,
  !!                      of same size as y_array.
  subroutine get_array_derivative(x_array, y_array, d_y_array_dx)
    real(dp), intent(in)       :: x_array(:), y_array(:)
    real(dp), intent(out)      :: d_y_array_dx(size(y_array))

    integer                    :: len_array, i

    len_array = size(y_array)

    do i = 2, len_array-1
      d_y_array_dx(i) = (y_array(i + 1) - y_array(i-1)) &
                      / 2.0d0*(x_array(i + 1) - x_array(i-1))
    end do

    ! TODO Figure out a way to set boundary conditions on the arrays,
    !      right now the values are copied into the boundary cells.
    ! left boundary condition
    d_y_array_dx(1) = d_y_array_dx(2)
    ! right boundary condition
    d_y_array_dx(len_array) = d_y_array_dx(len_array - 1)

  end subroutine get_array_derivative


  subroutine grid_clean()
    deallocate(grid)
    deallocate(grid_gauss)
  end subroutine grid_clean




end module mod_grid
