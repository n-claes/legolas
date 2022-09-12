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
  public :: get_accumulated_grid

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
    use mod_global_variables, only: geometry, x_start, x_end, &
      gauss_gridpts, dp_LIMIT, force_r0
    use mod_equilibrium_fields, only: initialise_equilibrium

    !> custom grid to use instead of the default one, optional
    real(dp), intent(in), optional  :: custom_grid(:)
    integer   :: i
    real(dp)  :: dx

    if (geometry == "") then
      call log_message("geometry must be set in submodule/parfile", level="error")
      return
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
      if (size(custom_grid) > gridpts) then
        call log_message( &
          "custom grid size larger than expected: "// str(size(custom_grid)) // &
          ". Grid will be trimmed down to " // str(gridpts) // " points", &
          level="info" &
        )
      else if (size(custom_grid) /= gridpts) then
        call log_message( &
          "custom grid: sizes do not match! Expected "// str(gridpts) // &
          " points but got " // str(size(custom_grid)), &
          level="error" &
        )
        return
      end if
      ! check if grid is monotonous
      do i = 1, size(custom_grid(1:gridpts)) - 1
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
      grid = custom_grid(1:gridpts)
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

    call set_grid_gauss()
    call set_scale_factor()
    call initialise_equilibrium()
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
        xi(j) = 0.5d0 * dx * gaussian_nodes(j) + 0.5d0 * (x_lo + x_hi)
        idx   = (i - 1) * n_gauss + j
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
      call log_message( &
        "geometry not defined correctly: " // trim(geometry), level="error" &
      )
    end if
  end subroutine set_scale_factor


  !> Generates a non-uniform grid based on a accumulation profile of choice.
  !! @note  The p1, p2, p3 and p4 parameters are not needed for all profiles.
  !!        If not needed, they can be set to NaN. See the individual
  !!        accumulation profiles to determine their necessity. @endnote
  subroutine get_accumulated_grid(accumulated_grid, grid_profile, p1, p2, p3, p4)
    use mod_global_variables, only: x_start, x_end, gridpts, set_gridpts, dp_LIMIT
    use mod_logging, only: log_message
    use mod_grid_profiles

    real(dp), intent(out) :: accumulated_grid(gridpts)
    character(len=*), intent(in) :: grid_profile
    real(dp), intent(in) :: p1, p2, p3, p4

    integer  :: i, new_gpts, idx(gridpts)
    real(dp) :: cx, dx, kappa, xbar(gridpts)

    abstract interface
      function func(x, p1, p2, p3, p4) result(dx)
        use mod_global_variables, only: dp
        real(dp) :: dx
        real(dp), intent(in) :: x, p1, p2, p3, p4
      end function func
    end interface

    ! pointer for the dx function, initialised to null
    procedure(func), pointer :: profile_dx => null()

    select case(grid_profile)
    case("gaussian")
      profile_dx => gaussian_dx
    case("sinusoidal")
      profile_dx => sine_dx
    case default
      call log_message( &
        "grid profile not recognised: " // trim(grid_profile), level="error" &
      )
    end select

    xbar(1) = x_start
    cx = x_start
    new_gpts = 1
    do while (cx < x_end)
      new_gpts = new_gpts + 1
      ! When defining new functions, also give them p1, p2, p3 and p4 arguments,
      ! even if they are not used.
      dx = profile_dx(cx, p1, p2, p3, p4)
      if (dx <= 0.0d0) then
        call log_message('dx in accumulated grid generation is negative !', level='error')
      end if
      cx = xbar(new_gpts-1) + dx
      xbar(new_gpts) = cx
    end do
    kappa = (x_end - xbar(new_gpts-1)) / (xbar(new_gpts) - xbar(new_gpts-1))

    accumulated_grid(1) = x_start
    do i = 1, new_gpts-1
      accumulated_grid(i+1) = xbar(i) + kappa * profile_dx(xbar(i), p1, p2, p3, p4)
    end do
    accumulated_grid(new_gpts) = x_end

    if (grid_profile == 'gaussian' .and. (x_start + x_end < dp_LIMIT) .and. abs(p2) < dp_LIMIT) then
      do i = 1, new_gpts
        if (accumulated_grid(i) > 0) then
          idx = i
          exit
        end if
      end do
      new_gpts = 2*idx(1)-1
      accumulated_grid(idx(1)) = 0.0d0
      do i = 1, idx(1)-1
        accumulated_grid(idx(1)+i) = -accumulated_grid(idx(1)-i)
      end do
    end if
    call set_gridpts(new_gpts)
  end subroutine get_accumulated_grid


  !> Cleanup routine, deallocates the arrays at module scope.
  subroutine grid_clean()
    deallocate(grid)
    deallocate(grid_gauss)
    deallocate(eps_grid)
    deallocate(d_eps_grid_dr)
  end subroutine grid_clean

end module mod_grid
