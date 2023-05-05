module mod_gravity_settings
  implicit none

  private

  type, public :: gravity_settings_t
    logical, private :: has_external_gravity

  contains

    procedure, public :: enable
    procedure, public :: disable
    procedure, public :: is_enabled
  end type gravity_settings_t

  public :: new_gravity_settings

contains

  pure function new_gravity_settings() result(gravity)
    type(gravity_settings_t) :: gravity
    gravity%has_external_gravity = .false.
  end function new_gravity_settings

  pure subroutine enable(this)
    class(gravity_settings_t), intent(inout) :: this
    this%has_external_gravity = .true.
  end subroutine enable

  pure subroutine disable(this)
    class(gravity_settings_t), intent(inout) :: this
    this%has_external_gravity = .false.
  end subroutine disable

  pure logical function is_enabled(this)
    class(gravity_settings_t), intent(in) :: this
    is_enabled = this%has_external_gravity
  end function is_enabled


end module mod_gravity_settings
