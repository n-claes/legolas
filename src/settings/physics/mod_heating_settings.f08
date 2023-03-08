module mod_heating_settings
  implicit none

  private

  type, public :: heating_settings_t
    logical, private :: has_heating

  contains

    procedure, public :: enable
    procedure, public :: disable
    procedure, public :: is_enabled
  end type heating_settings_t

  public :: new_heating_settings

contains

  pure function new_heating_settings() result(heating)
    type(heating_settings_t) :: heating
    heating%has_heating = .false.
  end function new_heating_settings


  pure logical function is_enabled(this)
    class(heating_settings_t), intent(in) :: this
    is_enabled = this%has_heating
  end function is_enabled


  pure subroutine enable(this)
    class(heating_settings_t), intent(inout) :: this
    this%has_heating = .true.
  end subroutine enable


  pure subroutine disable(this)
    class(heating_settings_t), intent(inout) :: this
    this%has_heating = .false.
  end subroutine disable

end module mod_heating_settings
