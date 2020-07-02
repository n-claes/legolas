! =============================================================================
!> @brief   Module containing all grid-related things.
!! @details Contains subroutines to create the base grid, Gaussian grid and
!!          scale factors. Also handles mesh accumulation.
!!          An integral of \e f(x) in <em>[a, b]</em> can be approximated with
!!          \f[ f(x) \approx 0.5(b-a)\sum_{i=1}^n\bigl[w_i f\bigl(0.5(b-a)x_i + 0.5(a+b)\bigr)\bigr] \f]
!!          where \f$w_i\f$ and \f$x_i\f$ are the weights and nodes of the Gaussian quadrature.
!!          The Gaussian grid is hence set up in every interval <em>[a, b]</em> across the
!!          four nodes \e j as
!!          \f[ x_i(j) = 0.5 * dx * w_i(j) + 0.5(a + b) \f]
module mod_grid
  use mod_global_variables, only: dp, gridpts
  use mod_logging, only: log_message
  implicit none

  private

  !> array containing the base grid
  real(dp), allocatable   :: grid(:)
  !> array containing the Gaussian grid
  real(dp), allocatable   :: grid_gauss(:)
  !> array containing the scale factor epsilon
  real(dp), allocatable   :: eps_grid(:)
  !> array containing the derivative of the scale factor epsilon
  real(dp), allocatable   :: d_eps_grid_dr(:)

  public :: grid
  public :: grid_gauss
  public :: eps_grid
  public :: d_eps_grid_dr

  public :: initialise_grid
  public :: set_grid_gauss
  public :: grid_clean

contains


  !> @brief   General grid initialisations.
  !! @details Initialises both the regular grid and the Gaussian grid.
  !!          Does calls to the mesh accumulation routines if needed,
  !!          and sets the scale factor and its derivative.
  !! @note  In order to avoid spurious eigenvalues in cylindrical geometry,
  !!        the default start is set at <tt> r = 0.025 </tt> instead of
  !!        <tt> r = 0 </tt>. The latter can be explicitly enforced by
  !!        setting \p force_r0 to \p True in the parfile.
  !! @warning   Explicitly forcing zero on-axis can give rise to spurious eigenvalues and
  !!            all kinds of other issues (huge values, division by zero, etc).
  !! @exception Error   If the geometry is not set.
  !! @exception Warning   If r=0 is forced in cylindrical geometry
  subroutine initialise_grid()
    use mod_global_variables, only: geometry, mesh_accumulation, x_start, x_end, &
                                    gauss_gridpts, dp_LIMIT, force_r0

    integer                  :: i
    real(dp)                 :: dx

    if (geometry == "") then
      call log_message("geometry must be set in submodule/parfile", level='error')
    end if

    allocate(grid(gridpts))
    allocate(grid_gauss(gauss_gridpts))
    allocate(eps_grid(gauss_gridpts))
    allocate(d_eps_grid_dr(gauss_gridpts))

    grid       = 0.0d0
    grid_gauss = 0.0d0

    if (geometry == 'cylindrical' .and. abs(x_start) < dp_LIMIT) then
      if (.not. force_r0) then
        x_start = 2.5d-2
      else
        call log_message("forcing on-axis r in cylindrical geometry. This may lead to spurious &
                        &real/imaginary eigenvalues.", level='warning')
        x_start = 0.0d0
      end if
    end if

    ! minus one here to include x_end
    dx = (x_end - x_start) / (gridpts-1)
    do i = 1, gridpts
      grid(i) = x_start + (i - 1)*dx
    end do

    if (mesh_accumulation) then
      call accumulate_mesh()
    end if

    call set_grid_gauss()
    call set_scale_factor()
  end subroutine initialise_grid


  !> Sets up grid_gauss, that is, the grid evaluated in the four
  !! Gaussian points.
  !> @brief   Sets up the Gaussian grid.
  !! @details Creates the Gaussian grid by evaluating the weights
  !!          at the four Gaussian nodes.
  subroutine set_grid_gauss()
    use mod_global_variables, only: gaussian_nodes, n_gauss

    real(dp)              :: x_lo, x_hi, dx, xi(n_gauss)
    integer               :: i, j, idx

    ! evaluates the nodes in the Gaussian quadrature
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


  !> @brief   Sets the scale factor and its derivative.
  !! @details The scale factor to switch between Cartesian and cylindrical geometries
  !!          is set here, along with its derivative. For cylindrical the scale factor
  !!          is simply equal to the Gaussian grid, and its derivative is unity.
  !!          For Cartesian the scale factor is unity and its derivative is zero.
  !! @exception Error If the geometry is not defined correctly.
  subroutine set_scale_factor()
    use mod_global_variables, only: geometry

    if (geometry == 'Cartesian') then
      eps_grid = 1.0d0
      d_eps_grid_dr = 0.0d0
    else if (geometry == 'cylindrical') then
      eps_grid = grid_gauss
      d_eps_grid_dr = 1.0d0
    else
      call log_message("geometry not defined correctly: " // trim(geometry), level='error')
    end if
  end subroutine set_scale_factor


  !> @brief   Does mesh accumulation.
  !! @details Re-grids the mesh to a non-uniform spacing using mesh accumulation.
  !!          This is based on two Gaussian curves with known widths, the grid is
  !!          accumulated near each Gaussian maximum using equidistribution based
  !!          on the integral under the curve defined by the function \p gaussian().
  subroutine accumulate_mesh()
    use mod_global_variables, only: x_start, x_end

    integer                  :: i, integral_gridpts
    integer                  :: integral_gridpts_1, integral_gridpts_2
    real(dp)                 :: dx, dx_0, xi, bgf, fact, dx_eq
    real(dp)                 :: gauss_xi, gauss_xi_eq
    real(dp)                 :: x_sum, x_sum_prev, x_norm
    real(dp)                 :: xi_weighted

    call log_message("redefining grid with mesh accumulation", level='info')

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



  !> @brief   Calculates a Gaussian curve based on a known width.
  !! @details Function to calculate a Gaussian curve based on known widths and
  !!          expected values. The Gaussian is evaluated in x.
  !! @param[in] x   value for x in the Gaussian function
  !! @param[in] bgf value for the background field, 0.3 by default
  !! @param[in] fact division factor in the Gaussian function, 1.0 by default
  !! @return f_gauss  value of the Gaussian function evaluated in x
  function gaussian(x, bgf, fact) result(f_gauss)
    use mod_global_variables, only: ev_1, ev_2, sigma_1, sigma_2
    use mod_physical_constants, only: dpi

    real(dp), intent(in)    :: x, bgf, fact
    real(dp)                :: f_gauss
    real(dp)                :: gauss_1, gauss_2, norm_1, norm_2

    norm_1 = 1.0d0 / (sigma_1 * sqrt(2*dpi))
    norm_2 = 1.0d0 / (sigma_2 * sqrt(2*dpi))

    gauss_1 = norm_1 * exp(-0.5d0 * ((x - ev_1) / sigma_1)**2)
    gauss_2 = norm_2 * exp(-0.5d0 * ((x - ev_2) / sigma_2)**2)

    f_gauss = bgf + (1.0d0 - bgf) * (gauss_1 + fact*gauss_2) / fact

  end function gaussian


  !> @brief   Cleanup routine
  !! @details Deallocates the arrays at module scope.
  subroutine grid_clean()
    deallocate(grid)
    deallocate(grid_gauss)
    deallocate(eps_grid)
    deallocate(d_eps_grid_dr)
  end subroutine grid_clean

end module mod_grid
