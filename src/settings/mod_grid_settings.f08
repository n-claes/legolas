module mod_grid_settings
  use mod_global_variables, only: dp, n_gauss
  implicit none

  private

  type, public :: grid_settings_t
    character(:), private, allocatable :: geometry
    integer, private :: gridpts
    integer, private :: gauss_gridpts
    integer, private :: ef_gridpts
    real(dp), private :: grid_start
    real(dp), private :: grid_end
    logical :: coaxial
    logical :: force_r0

  contains

    procedure, public :: set_geometry
    procedure, public :: get_geometry
    procedure, public :: set_gridpts
    procedure, public :: get_gridpts
    procedure, public :: get_gauss_gridpts
    procedure, public :: get_ef_gridpts
    procedure, public :: set_grid_boundaries
    procedure, public :: get_grid_start
    procedure, public :: get_grid_end
    procedure, public :: delete
  end type grid_settings_t

  public :: new_grid_settings

contains

  pure function new_grid_settings() result(grid_settings)
    type(grid_settings_t) :: grid_settings

    call grid_settings%set_geometry("Cartesian")
    call grid_settings%set_gridpts(50)
    call grid_settings%set_grid_boundaries(0.0_dp, 1.0_dp)
    grid_settings%coaxial = .false.
    grid_settings%force_r0 = .false.
  end function new_grid_settings


  pure subroutine set_geometry(this, geometry)
    class(grid_settings_t), intent(inout) :: this
    character(len=*), intent(in) :: geometry
    this%geometry = geometry
  end subroutine set_geometry


  pure function get_geometry(this) result(geometry)
    class(grid_settings_t), intent(in) :: this
    character(len=:), allocatable :: geometry
    geometry = trim(adjustl(this%geometry))
  end function get_geometry


  pure subroutine set_gridpts(this, gridpts)
    class(grid_settings_t), intent(inout) :: this
    integer, intent(in) :: gridpts
    this%gridpts = gridpts
    this%gauss_gridpts = n_gauss * (gridpts - 1)
    this%ef_gridpts = 2 * gridpts - 1
  end subroutine set_gridpts


  pure integer function get_gridpts(this)
    class(grid_settings_t), intent(in) :: this
    get_gridpts = this%gridpts
  end function get_gridpts


  pure integer function get_gauss_gridpts(this)
    class(grid_settings_t), intent(in) :: this
    get_gauss_gridpts = this%gauss_gridpts
  end function get_gauss_gridpts


  pure integer function get_ef_gridpts(this)
    class(grid_settings_t), intent(in) :: this
    get_ef_gridpts = this%ef_gridpts
  end function get_ef_gridpts


  pure subroutine set_grid_boundaries(this, grid_start, grid_end)
    class(grid_settings_t), intent(inout) :: this
    real(dp), intent(in) :: grid_start
    real(dp), intent(in) :: grid_end

    this%grid_start = grid_start
    this%grid_end = grid_end
  end subroutine set_grid_boundaries


  pure real(dp) function get_grid_start(this)
    class(grid_settings_t), intent(in) :: this
    get_grid_start = this%grid_start
  end function get_grid_start


  pure real(dp) function get_grid_end(this)
    class(grid_settings_t), intent(in) :: this
    get_grid_end = this%grid_end
  end function get_grid_end


  pure subroutine delete(this)
    class(grid_settings_t), intent(inout) :: this
    if (allocated(this%geometry)) deallocate(this%geometry)
  end subroutine delete

end module mod_grid_settings
