module mod_settings
  use mod_global_variables, only: str_len_arr
  use mod_dims, only: dims_t, new_block_dims
  use mod_io_settings, only: io_t, new_io_settings
  use mod_solver_settings, only: solvers_t, new_solver_settings
  implicit none

  private

  type, public :: settings_t
    character(len=:), private, allocatable :: state_vector(:)
    character(len=:), private, allocatable :: physics_type
    type(dims_t), public :: dims
    type(io_t), public :: io
    type(solvers_t), public :: solvers

  contains

    procedure, public :: initialise
    procedure, public :: get_state_vector
    procedure, public :: get_physics_type
    procedure, public :: delete
  end type settings_t

  public :: new_settings

contains

  pure function new_settings() result(settings)
    type(settings_t) :: settings
  end function new_settings


  pure subroutine initialise(this, physics_type, gridpts)
    class(settings_t), intent(inout) :: this
    character(len=*), intent(in) :: physics_type
    integer, intent(in) :: gridpts

    this%physics_type = physics_type
    call set_state_vector(this, physics_type)
    call this%dims%set_block_dims(nb_eqs=size(this%state_vector), gridpts=gridpts)
    this%io = new_io_settings()
    this%solvers = new_solver_settings()
  end subroutine initialise


  pure subroutine set_state_vector(settings, physics_type)
    type(settings_t), intent(inout) :: settings
    character(len=*), intent(in) :: physics_type

    select case(physics_type)
      case("hd")
        settings%state_vector = [ &
          character(len=str_len_arr) :: "rho", "v1", "v2", "v3", "T" &
        ]
      case default
        settings%state_vector = [ &
          character(len=str_len_arr) :: "rho", "v1", "v2", "v3", "T", "a1", "a2", "a3" &
        ]
      end select
  end subroutine set_state_vector


  pure function get_state_vector(this) result(state_vector)
    class(settings_t), intent(in) :: this
    character(len=:), allocatable :: state_vector(:)

    state_vector = this%state_vector
  end function get_state_vector


  pure function get_physics_type(this) result(physics_type)
    class(settings_t), intent(in) :: this
    character(len=:), allocatable :: physics_type

    physics_type = this%physics_type
  end function get_physics_type


  pure subroutine delete(this)
    class(settings_t), intent(inout) :: this

    if (allocated(this%state_vector)) deallocate(this%state_vector)
    if (allocated(this%physics_type)) deallocate(this%physics_type)
    call this%io%delete()
    call this%solvers%delete()
  end subroutine delete

end module mod_settings
