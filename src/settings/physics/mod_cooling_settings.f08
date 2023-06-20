module mod_cooling_settings
  implicit none

  private

  type, public :: cooling_settings_t
    integer, private :: n_interp
    character(:), private, allocatable :: cooling_curve
    logical, private :: has_cooling

  contains

    procedure, public :: enable
    procedure, public :: disable
    procedure, public :: is_enabled
    procedure, public :: set_interpolation_points
    procedure, public :: get_interpolation_points
    procedure, public :: set_cooling_curve
    procedure, public :: get_cooling_curve
  end type cooling_settings_t

  public :: new_cooling_settings

contains

  pure function new_cooling_settings() result(cooling)
    use mod_cooling_curve_names, only: NOTHING
    type(cooling_settings_t) :: cooling

    cooling%has_cooling = .false.
    call cooling%set_cooling_curve(NOTHING)
    call cooling%set_interpolation_points(4000)
  end function new_cooling_settings


  pure logical function is_enabled(this)
    class(cooling_settings_t), intent(in) :: this
    is_enabled = this%has_cooling
  end function is_enabled


  pure subroutine enable(this)
    class(cooling_settings_t), intent(inout) :: this
    this%has_cooling = .true.
  end subroutine enable


  pure subroutine disable(this)
    class(cooling_settings_t), intent(inout) :: this
    this%has_cooling = .false.
  end subroutine disable


  pure subroutine set_interpolation_points(this, n_interp)
    class(cooling_settings_t), intent(inout) :: this
    integer, intent(in) :: n_interp
    this%n_interp = n_interp
  end subroutine set_interpolation_points


  pure integer function get_interpolation_points(this)
    class(cooling_settings_t), intent(in) :: this
    get_interpolation_points = this%n_interp
  end function get_interpolation_points


  pure subroutine set_cooling_curve(this, cooling_curve)
    class(cooling_settings_t), intent(inout) :: this
    character(len=*), intent(in) :: cooling_curve
    this%cooling_curve = cooling_curve
  end subroutine set_cooling_curve


  pure function get_cooling_curve(this) result(cooling_curve)
    class(cooling_settings_t), intent(in) :: this
    character(:), allocatable :: cooling_curve
    cooling_curve = trim(adjustl(this%cooling_curve))
  end function get_cooling_curve

end module mod_cooling_settings
