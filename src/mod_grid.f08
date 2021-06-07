! =============================================================================
!> Module containing all grid-related things.
!! Contains subroutines to create the base grid, Gaussian grid and
!! scale factors. Also handles mesh accumulation.
!! An integral of \(f(x)\) in \([a, b]\) can be approximated with
!! $$ f(x) \approx 0.5(b-a)\sum_{i=1}^n\bigl[w_i f\bigl(0.5(b-a)x_i + 0.5(a+b)\bigr)\bigr] $$
!! where \(w_i\) and \(x_i\) are the weights and nodes of the Gaussian quadrature.
!! The Gaussian grid is hence set up in every interval \([a, b]\) across the
!! four nodes \(j\) as
!! $$ x_i(j) = 0.5 * dx * w_i(j) + 0.5(a + b) $$
module mod_grid
  use mod_global_variables, only: dp, gridpts
  use mod_logging, only: log_message, str
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


  !> General grid initialisations.
  !! Initialises both the regular grid and the Gaussian grid.
  !! Does calls to the mesh accumulation routines if needed,
  !! and sets the scale factor and its derivative.
  !! @note  In order to avoid spurious eigenvalues in cylindrical geometry,
  !!        the default start is set at <tt> r = 0.025 </tt> instead of
  !!        <tt> r = 0 </tt>. The latter can be explicitly enforced by
  !!        setting <tt>force_r0</tt> to <tt>True</tt> in the parfile. @endnote
  !! @warning   Explicitly forcing zero on-axis can give rise to spurious eigenvalues and
  !!            all kinds of other issues (huge values, division by zero, etc). @endwarning
  !! @warning   Throws an error if the geometry is not set. @endwarning
  !! @warning   Throws a warning if <tt>r = 0</tt> is forced in cylindrical geometry.
  subroutine initialise_grid(custom_grid)
    use mod_global_variables, only: geometry, mesh_accumulation, x_start, x_end, &
      gauss_gridpts, dp_LIMIT, force_r0, equilibrium_type

    !> custom grid to use instead of the default one, optional
    real(dp), intent(in), optional  :: custom_grid(:)
    integer   :: i
    real(dp)  :: dx

    if (geometry == "") then
      call log_message("geometry must be set in submodule/parfile", level="error")
    end if

    allocate(grid(gridpts))
    allocate(grid_gauss(gauss_gridpts))
    allocate(eps_grid(gauss_gridpts))
    allocate(d_eps_grid_dr(gauss_gridpts))

    grid = 0.0d0
    grid_gauss = 0.0d0

    if (present(custom_grid)) then
      call log_message("using a custom base grid", level="info")
      ! check if size matches gridpoints
      if (size(custom_grid) /= gridpts) then
        call log_message( &
          "custom grid: sizes do not match! Expected "// str(gridpts) // &
          " points but got " // str(size(custom_grid)), &
          level="error" &
        )
        return
      end if
      ! check if grid is monotonous
      do i = 1, size(custom_grid) - 1
        if (.not. (custom_grid(i + 1) > custom_grid(i))) then
          call log_message( &
            "custom grid: supplied array is not monotone! Got x=" // &
            str(custom_grid(i)) // " at index " // str(i) // " and x=" // &
            str(custom_grid(i + 1)) // " at index " // str(i + 1), &
            level="error" &
          )
          return
        end if
      end do
      grid = custom_grid
    else
      if (geometry == "cylindrical" .and. abs(x_start) < dp_LIMIT) then
        if (.not. force_r0) then
          x_start = 2.5d-2
        else
          call log_message( &
            "forcing on-axis r in cylindrical geometry. This may lead to spurious &
            &real/imaginary eigenvalues.", &
            level="warning" &
          )
          x_start = 0.0d0
        end if
      end if
      
      ! minus one here to include x_end
      dx = (x_end - x_start) / (gridpts-1)
      do i = 1, gridpts
        grid(i) = x_start + (i - 1)*dx
      end do
    end if

    if (mesh_accumulation) then
      if (equilibrium_type == 'photospheric_flux_tube' .or. equilibrium_type == 'coronal_flux_tube') then
        call accumulate_mesh_special_fluxtube()
      else
        call accumulate_mesh()
      end if
    end if

    call set_grid_gauss()
    call set_scale_factor()
  end subroutine initialise_grid


  !> Sets up grid_gauss, that is, the grid evaluated in the four
  !! Gaussian points. This is done by evaluating the weights
  !! at the four Gaussian nodes.
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


  !> The scale factor to switch between Cartesian and cylindrical geometries
  !! is set here, along with its derivative. For cylindrical the scale factor
  !! is simply equal to the Gaussian grid, and its derivative is unity.
  !! For Cartesian the scale factor is unity and its derivative is zero.
  !! @warning   Throws an error if the geometry is not defined correctly.
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


  !> Special mesh accumulation routine for flux tubes. Instead of a Gaussian
  !! a 60-30-10 percentage-distribution is used, meaning that 60% of the gridpoints
  !! are located between the axis and \((a - width/2)\), 30% of the gridpoints are
  !! located at a \(width=0.1\) centered at the inner fluxtube radius \(a\),
  !! and the remaining 10% is located in the outside region.
  subroutine accumulate_mesh_special_fluxtube()
    use mod_global_variables, only: x_start, x_end
    use mod_equilibrium_params, only: r0

    real(dp)  :: width, a_l, a_r
    real(dp)  :: pct1, pct2, pct3
    real(dp)  :: dx, dx1, dx2, dx3
    integer   :: i, N1, N2, N3

    ! reset the grid
    grid = 0.0d0
    ! width of transition region
    width = 0.1d0
    ! start of transition region
    a_l = r0 - width / 2.0d0
    ! end of transition region
    a_r = a_l + width

    pct1 = 0.6d0  ! percentage of points in inner tube
    pct2 = 0.3d0  ! percentage of points in transition region
    pct3 = 0.1d0  ! percentage of points in outer tube

    N1 = int(pct1 * gridpts)
    dx1 = (a_l - x_start) / (N1 - 1)  ! -1 since first point is x0
    N2 = int(pct2 * gridpts)
    dx2 = (a_r - a_l) / N2
    N3 = int(pct3 * gridpts)
    dx3 = (x_end - a_r) / N3

    ! fill the grid
    grid(1) = x_start
    do i = 2, gridpts
      if (i <= N1) then
        dx = dx1
      else if (N1 < i .and. i <= N1 + N2) then
        dx = dx2
      else
        dx = dx3
      end if
      grid(i) = grid(i - 1) + dx
    end do
    call log_message("special mesh refinement for flux tubes, 60%-30%-10%", level='info')
  end subroutine accumulate_mesh_special_fluxtube


  !> Re-grids the mesh to a non-uniform spacing using mesh accumulation.
  !! This is based on two Gaussian curves with known widths, the grid is
  !! accumulated near each Gaussian maximum using equidistribution based
  !! on the integral under the curve defined by the function <tt>gaussian()</tt>.
  !! @todo  This method will most likely change in the future.
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



  !> Function to calculate a Gaussian curve based on known widths and
  !! expected values. Returns the value of the Gaussian function,
  !! evaluated in x.
  !! @todo  This method will most likely change in the future.
  function gaussian(x, bgf, fact) result(f_gauss)
    use mod_global_variables, only: ev_1, ev_2, sigma_1, sigma_2
    use mod_physical_constants, only: dpi

    !> value for x in the Gaussian function
    real(dp), intent(in)    :: x
    !> value for the background field, 0.3 by default
    real(dp), intent(in)    :: bgf
    !> division factor in the Gaussian function, 1.0 by default
    real(dp), intent(in)    :: fact
    real(dp)                :: f_gauss
    real(dp)                :: gauss_1, gauss_2, norm_1, norm_2

    norm_1 = 1.0d0 / (sigma_1 * sqrt(2*dpi))
    norm_2 = 1.0d0 / (sigma_2 * sqrt(2*dpi))

    gauss_1 = norm_1 * exp(-0.5d0 * ((x - ev_1) / sigma_1)**2)
    gauss_2 = norm_2 * exp(-0.5d0 * ((x - ev_2) / sigma_2)**2)

    f_gauss = bgf + (1.0d0 - bgf) * (gauss_1 + fact*gauss_2) / fact

  end function gaussian


  !> Cleanup routine, deallocates the arrays at module scope.
  subroutine grid_clean()
    deallocate(grid)
    deallocate(grid_gauss)
    deallocate(eps_grid)
    deallocate(d_eps_grid_dr)
  end subroutine grid_clean

end module mod_grid
