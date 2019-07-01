!
! MODULE: mod_grid
!
!> @author
!> Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Module to create the grid and to do mesh accumulation if needed.
!
module mod_grid
  use mod_global_variables
  implicit none

  public

  !> Array containing the regular (coarse) grid
  real(dp), allocatable      :: grid(:)
  !> New array with 4x the length of grid, due to the 4 nodes of the
  !! Gaussian quadrature
  real(dp), allocatable      :: grid_gauss(:)

  private :: accumulate_mesh
  private :: gaussian

contains

  !> Initialises both the coarse grid and grid_gauss. Does mesh accumulation
  !! if desired.
  subroutine initialise_grid()
    integer                  :: i
    real(dp)                 :: dx

    allocate(grid(gridpts))
    allocate(grid_gauss(gauss_gridpts))

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

    call set_grid_gauss()

  end subroutine initialise_grid

  !> Sets up grid_gauss, that is, the grid evaluated in the four
  !! Gaussian points.
  subroutine set_grid_gauss()
    real(dp)              :: x_lo, x_hi, dx, xi(n_gauss)
    integer               :: i, j, idx

    ! Check for origin in cylindrical coordinates
    ! \TODO: is this needed? This is automatically handled in the loop below??
    if (geometry == "cylindrical" .and. grid_gauss(1) .le. 1.0d-5) then
      grid_gauss(1) = 1.0d-5
    end if

    ! Evaluate grid_gauss in nodes of Gaussian quadrature.
    ! An integral of f(x) in [a, b] can be approximated by
    ! 0.5*(b-a) * SUM[i from 1 -> n] ( wi * f( 0.5*(b-a)*xi + 0.5*(b-a)) )
    ! where wi and xi are the weights and nodes at i (so 1 to 4 here).
    ! Hence we need the gridpoints equal to the evaluation points of f(x).
    do i = 1, gridpts - 1
      x_lo = grid(i)
      x_hi = grid(i + 1)
      dx   = x_hi - x_lo

      do j = 1, n_gauss
        xi(j) = 0.5 * dx * gaussian_nodes(j) + 0.5*(x_lo + x_hi)
        idx   = (i - 1)*n_gauss + j
        grid_gauss(idx) = xi(j)
      end do
    end do
  end subroutine set_grid_gauss

  !> Subroutine to re-grid the mesh to a non-uniform spacing.
  !! This is based on two Gaussian curves with known widths (from
  !! mod_global_variables); the grid is accumulated near each Gaussian maximum,
  !! using equidistribution based on the integral under the curve defined
  !! by the function gaussian().
  subroutine accumulate_mesh()
    integer                  :: i, integral_gridpts
    integer                  :: integral_gridpts_1, integral_gridpts_2
    real(dp)                 :: dx, dx_0, xi, bgf, fact, dx_eq
    real(dp)                 :: gauss_xi, gauss_xi_eq
    real(dp)                 :: x_sum, x_sum_prev, x_norm
    real(dp)                 :: xi_weighted

    print*,"Redefining grid with mesh accumulation"

    integral_gridpts = gridpts - 1

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
  !! @param[in] x   Value for x in the Gaussian function
  !! @param[in] bgf Value for the background field, hardcoded to 0.3
  !! @param[in] fact Division factor in the Gaussian function, hardcoded to 1.0
  !! @return f_gauss  Value of the Gaussian function, evaluated in x
  function gaussian(x, bgf, fact) result(f_gauss)
    use mod_physical_constants

    real(dp), intent(in)    :: x, bgf, fact
    real(dp)                :: f_gauss
    real(dp)                :: gauss_1, gauss_2, norm_1, norm_2

    norm_1 = 1.0d0 / (sigma_1 * sqrt(2*dpi))
    norm_2 = 1.0d0 / (sigma_2 * sqrt(2*dpi))

    gauss_1 = norm_1 * exp(-0.5d0 * ((x - ev_1) / sigma_1)**2)
    gauss_2 = norm_2 * exp(-0.5d0 * ((x - ev_2) / sigma_2)**2)

    f_gauss = bgf + (1.0d0 - bgf) * (gauss_1 + fact*gauss_2) / fact

  end function gaussian

  !> Deallocates arrays defined in this module.
  subroutine grid_clean()
    deallocate(grid)
    deallocate(grid_gauss)
  end subroutine grid_clean




end module mod_grid
