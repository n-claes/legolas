! =============================================================================
!> Module containing all grid-related things.
!! Contains subroutines to create the base grid, Gaussian grid and
!! scale factors.
!! An integral of \(f(x)\) in \([a, b]\) can be approximated with
!! $$
!! f(x) \approx 0.5(b-a)\sum_{i=1}^n\bigl[w_i f\bigl(0.5(b-a)x_i + 0.5(a+b)\bigr)\bigr] !! $$
!! where \(w_i\) and \(x_i\) are the weights and nodes of the Gaussian quadrature.
!! The Gaussian grid is hence set up in every interval \([a, b]\) across the
!! nodes \(j\) as
!! $$ x_i(j) = 0.5 * dx * w_i(j) + 0.5(a + b) $$
module mod_grid
  use mod_global_variables, only: dp, NaN
  use mod_logging, only: logger, str
  use mod_settings, only: settings_t
  implicit none

  private

  type, public :: grid_t
    real(dp), allocatable, public :: base_grid(:)
    real(dp), allocatable, public :: gaussian_grid(:)
    real(dp), allocatable, public :: ef_grid(:)

    type(settings_t), pointer, private :: settings
    logical, private :: is_initialised
    logical, private :: uses_custom_base_grid

  contains

    procedure, private :: set_base_grid
    procedure, private :: set_gaussian_grid
    procedure, private :: set_ef_grid

    procedure, public :: initialise
    procedure, public :: set_custom_grid
    procedure, public :: get_eps
    procedure, public :: get_deps
    procedure, public :: delete
  end type grid_t

  public :: new_grid

contains

  function new_grid(settings) result(grid)
    type(settings_t), target, intent(in) :: settings
    type(grid_t) :: grid
    grid%settings => settings
    grid%is_initialised = .false.
    grid%uses_custom_base_grid = .false.
  end function new_grid


  subroutine initialise(this)
    class(grid_t), intent(inout) :: this

    if (this%is_initialised) return

    if (.not. this%uses_custom_base_grid) call this%set_base_grid()
    call this%set_gaussian_grid()
    if (this%settings%io%write_eigenfunctions) call this%set_ef_grid()
    this%is_initialised = .true.
  end subroutine initialise


  subroutine set_custom_grid(this, custom)
    class(grid_t), intent(inout) :: this
    real(dp), intent(in) :: custom(:)

    if (.not. is_valid_custom_base_grid(this%settings, custom)) return
    allocate(this%base_grid, source=custom)
    this%uses_custom_base_grid = .true.
  end subroutine set_custom_grid


  pure subroutine set_base_grid(this)
    class(grid_t), intent(inout) :: this
    integer :: i, gridpts
    real(dp) :: dx, grid_start, grid_end

    gridpts = this%settings%grid%get_gridpts()
    allocate(this%base_grid(gridpts))
    grid_start = this%settings%grid%get_grid_start()
    grid_end = this%settings%grid%get_grid_end()
    ! minus one to include the end point
    dx = (grid_end - grid_start) / (gridpts - 1)
    do i = 1, gridpts
      this%base_grid(i) = grid_start + (i - 1) * dx
    end do
  end subroutine set_base_grid


  pure subroutine set_gaussian_grid(this)
    use mod_global_variables, only: gaussian_nodes, n_gauss

    class(grid_t), intent(inout) :: this
    integer :: i, j, gauss_idx
    real(dp) :: x_lo, x_hi, dx

    allocate(this%gaussian_grid(this%settings%grid%get_gauss_gridpts()))
    do i = 1, size(this%base_grid) - 1
      x_lo = this%base_grid(i)
      x_hi = this%base_grid(i + 1)
      dx = x_hi - x_lo
      do j = 1, n_gauss
        gauss_idx = (i - 1) * n_gauss + j
        this%gaussian_grid(gauss_idx) = ( &
          0.5_dp * dx * gaussian_nodes(j) + 0.5_dp * (x_lo + x_hi) &
        )
      end do
    end do
  end subroutine set_gaussian_grid


  pure subroutine set_ef_grid(this)
    class(grid_t), intent(inout) :: this
    integer :: idx

    allocate(this%ef_grid(this%settings%grid%get_ef_gridpts()))
    ! first gridpoint, left edge
    this%ef_grid(1) = this%base_grid(1)
    ! other gridpoints
    do idx = 1, this%settings%grid%get_gridpts() - 1
      ! position of center point in grid interval
      this%ef_grid(2 * idx) = 0.5_dp * (this%base_grid(idx) + this%base_grid(idx + 1))
      ! position of end point in grid interval
      this%ef_grid(2 * idx + 1) = this%base_grid(idx + 1)
    end do
  end subroutine set_ef_grid


  impure real(dp) elemental function get_eps(this, x)
    class(grid_t), intent(in) :: this
    real(dp), intent(in) :: x

    select case(this%settings%grid%get_geometry())
    case("Cartesian")
      get_eps = 1.0_dp
    case("cylindrical")
      get_eps = x
    case default
      get_eps = NaN
      call logger%error("geometry has no defined scale factor")
    end select
  end function get_eps


  impure real(dp) elemental function get_deps(this)
    class(grid_t), intent(in) :: this

    select case(this%settings%grid%get_geometry())
    case("Cartesian")
      get_deps = 0.0_dp
    case("cylindrical")
      get_deps = 1.0_dp
    case default
      get_deps = NaN
      call logger%error("geometry has no defined scale factor derivative")
    end select
  end function get_deps


  pure subroutine delete(this)
    class(grid_t), intent(inout) :: this
    if (allocated(this%base_grid)) deallocate(this%base_grid)
    if (allocated(this%gaussian_grid)) deallocate(this%gaussian_grid)
    if (allocated(this%ef_grid)) deallocate(this%ef_grid)
    nullify(this%settings)
    this%is_initialised = .false.
  end subroutine delete


  logical function is_valid_custom_base_grid(settings, custom_grid)
    type(settings_t), intent(inout) :: settings
    real(dp), intent(in) :: custom_grid(:)
    integer :: i, gridpts

    is_valid_custom_base_grid = .false.
    gridpts = settings%grid%get_gridpts()
    if (size(custom_grid) /= gridpts) then
      call logger%error( &
        "custom grid: sizes do not match! Expected "// str(gridpts) // &
        " points but got " // str(size(custom_grid)) &
      )
      return
    end if
    ! check monotonicity
    do i = 1, size(custom_grid) - 1
      if (.not. (custom_grid(i + 1) > custom_grid(i))) then
        call logger%error( &
          "custom grid: supplied array is not monotone! Got x=" // &
          str(custom_grid(i)) // " at index " // str(i) // " and x=" // &
          str(custom_grid(i + 1)) // " at index " // str(i + 1) &
        )
        return
      end if
    end do
    ! ensure grid start/end are consistent
    call settings%grid%set_grid_boundaries(custom_grid(1), custom_grid(gridpts))
    is_valid_custom_base_grid = .true.
  end function is_valid_custom_base_grid

end module mod_grid
