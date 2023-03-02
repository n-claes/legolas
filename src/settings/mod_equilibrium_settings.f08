module mod_equilibrium_settings
  implicit none

  type, public :: equilibrium_settings_t
    character(:), allocatable :: equilibrium_type
    character(:), allocatable :: boundary_type
    logical :: use_defaults

  contains

    procedure, public :: set_equilibrium_type
    procedure, public :: get_equilibrium_type
    procedure, public :: set_boundary_type
    procedure, public :: get_boundary_type
    procedure, public :: delete
  end type equilibrium_settings_t

  public :: new_equilibrium_settings

contains

  pure function new_equilibrium_settings() result(equilibrium)
    type(equilibrium_settings_t) :: equilibrium

    call equilibrium%set_equilibrium_type("adiabatic_homo")
    call equilibrium%set_boundary_type("wall")
    equilibrium%use_defaults = .true.
  end function new_equilibrium_settings


  pure subroutine set_equilibrium_type(this, equilibrium_type)
    class(equilibrium_settings_t), intent(inout) :: this
    character(len=*), intent(in) :: equilibrium_type
    this%equilibrium_type = equilibrium_type
  end subroutine set_equilibrium_type


  pure function get_equilibrium_type(this) result(equilibrium_type)
    class(equilibrium_settings_t), intent(in) :: this
    character(len=:), allocatable :: equilibrium_type
    equilibrium_type = this%equilibrium_type
  end function get_equilibrium_type


  pure subroutine set_boundary_type(this, boundary_type)
    class(equilibrium_settings_t), intent(inout) :: this
    character(len=*), intent(in) :: boundary_type
    this%boundary_type = boundary_type
  end subroutine set_boundary_type


  pure function get_boundary_type(this) result(boundary_type)
    class(equilibrium_settings_t), intent(in) :: this
    character(len=:), allocatable :: boundary_type
    boundary_type = this%boundary_type
  end function get_boundary_type

  pure subroutine delete(this)
    class(equilibrium_settings_t), intent(inout) :: this
    if (allocated(this%equilibrium_type)) deallocate(this%equilibrium_type)
    if (allocated(this%boundary_type)) deallocate(this%boundary_type)
  end subroutine delete

end module mod_equilibrium_settings
