module mod_flow_settings
  implicit none

  type, public :: flow_settings_t
    logical, private :: has_flow

  contains

    procedure, public :: enable
    procedure, public :: disable
    procedure, public :: is_enabled
  end type flow_settings_t

  public :: new_flow_settings

contains

  pure function new_flow_settings() result(flow)
    type(flow_settings_t) :: flow
    flow%has_flow = .false.
  end function new_flow_settings


  pure logical function is_enabled(this)
    class(flow_settings_t), intent(in) :: this
    is_enabled = this%has_flow
  end function is_enabled


  pure subroutine enable(this)
    class(flow_settings_t), intent(inout) :: this
    this%has_flow = .true.
  end subroutine enable


  pure subroutine disable(this)
    class(flow_settings_t), intent(inout) :: this
    this%has_flow = .false.
  end subroutine disable

end module mod_flow_settings
