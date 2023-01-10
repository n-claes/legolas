module mod_hall_settings
  use mod_global_variables, only: dp
  implicit none

  private

  type, public :: hall_settings_t
    logical, private :: has_hall
    logical, private :: use_hall_substitution
    logical, private :: electron_inertia
    real(dp), private :: electron_fraction

    logical, public :: use_dropoff
    logical, public :: use_inertia_dropoff

  contains

    procedure, public :: enable
    procedure, public :: enable_electron_inertia
    procedure, public :: disable
    procedure, public :: is_enabled
    procedure, public :: is_using_substitution
    procedure, public :: has_electron_inertia
    procedure, public :: set_electron_fraction
    procedure, public :: get_electron_fraction
  end type hall_settings_t

  public :: new_hall_settings

contains

  pure function new_hall_settings() result(hall)
    type(hall_settings_t) :: hall

    hall%has_hall = .false.
    hall%use_hall_substitution = .true.
    hall%electron_inertia = .false.
    hall%electron_fraction = 0.5_dp
    hall%use_dropoff = .false.
    hall%use_inertia_dropoff = .false.
  end function new_hall_settings


  pure logical function is_enabled(this)
    class(hall_settings_t), intent(in) :: this
    is_enabled = this%has_hall
  end function is_enabled


  pure logical function is_using_substitution(this)
    class(hall_settings_t), intent(in) :: this
    is_using_substitution = this%use_hall_substitution
  end function is_using_substitution


  pure subroutine enable(this)
    class(hall_settings_t), intent(inout) :: this
    this%has_hall = .true.
    ! we always use substitution but keep this for future use
    this%use_hall_substitution = .true.
  end subroutine enable


  pure subroutine enable_electron_inertia(this)
    class(hall_settings_t), intent(inout) :: this
    this%electron_inertia = .true.
  end subroutine enable_electron_inertia


  pure subroutine disable(this)
    class(hall_settings_t), intent(inout) :: this
    this%has_hall = .false.
    this%electron_inertia = .false.
  end subroutine disable


  pure logical function has_electron_inertia(this)
    class(hall_settings_t), intent(in) :: this
    has_electron_inertia = this%electron_inertia
  end function has_electron_inertia


  pure subroutine set_electron_fraction(this, electron_fraction)
    class(hall_settings_t), intent(inout) :: this
    real(dp), intent(in) :: electron_fraction
    this%electron_fraction = electron_fraction
  end subroutine set_electron_fraction


  pure real(dp) function get_electron_fraction(this)
    class(hall_settings_t), intent(in) :: this
    get_electron_fraction = this%electron_fraction
  end function get_electron_fraction

end module mod_hall_settings
