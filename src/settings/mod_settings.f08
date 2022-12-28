module mod_settings
  use mod_global_variables, only: str_len_arr
  use mod_dims, only: dims_t, new_block_dims
  use mod_io_settings, only: io_t, new_io_settings
  use mod_solver_settings, only: solvers_t, new_solver_settings
  use mod_physics_settings, only: physics_t, new_physics_settings
  use mod_grid_settings, only: grid_settings_t, new_grid_settings
  use mod_equilibrium_settings, only: equilibrium_settings_t, new_equilibrium_settings
  use mod_units, only: units_t, new_unit_system
  implicit none

  private

  type, public :: settings_t
    ! note: weird gfortran 8 bug here when using (len=:) for state_vector.
    ! This sometimes leads to wrong array allocation where every entry equals the
    ! one at the last index? Unable to reproduce with compiler versions >8.
    character(len=str_len_arr), private, allocatable :: state_vector(:)
    character(len=:), private, allocatable :: physics_type
    integer, private :: nb_eqs
    type(dims_t), public :: dims
    type(io_t), public :: io
    type(solvers_t), public :: solvers
    type(physics_t), public :: physics
    type(grid_settings_t), public :: grid
    type(equilibrium_settings_t), public :: equilibrium
    type(units_t), public :: units

  contains

    procedure, public :: set_state_vector
    procedure, public :: get_state_vector
    procedure, public :: state_vector_is_set
    procedure, public :: get_physics_type
    procedure, public :: get_nb_eqs
    procedure, public :: update_block_dimensions
    procedure, public :: delete

    procedure, private :: set_nb_eqs
  end type settings_t

  public :: new_settings

contains

  pure function new_settings() result(settings)
    type(settings_t) :: settings

    call settings%set_state_vector(physics_type="mhd")
    settings%dims = new_block_dims()
    settings%io = new_io_settings()
    settings%solvers = new_solver_settings()
    settings%physics = new_physics_settings()
    settings%grid = new_grid_settings()
    settings%equilibrium = new_equilibrium_settings()
    settings%units = new_unit_system()
    call settings%update_block_dimensions()
  end function new_settings


  pure subroutine set_state_vector(this, physics_type)
    class(settings_t), intent(inout) :: this
    character(len=*), intent(in) :: physics_type

    this%physics_type = physics_type
    select case(physics_type)
      case("hd")
        this%state_vector = [character(3) :: "rho", "v1", "v2", "v3", "T"]
      case default
        this%state_vector = [ &
          character(3) :: "rho", "v1", "v2", "v3", "T", "a1", "a2", "a3" &
        ]
      end select
    call this%set_nb_eqs(size(this%state_vector))
  end subroutine set_state_vector


  pure logical function state_vector_is_set(this)
    class(settings_t), intent(in) :: this

    state_vector_is_set = allocated(this%state_vector)
  end function state_vector_is_set


  pure function get_state_vector(this) result(state_vector)
    class(settings_t), intent(in) :: this
    character(len=:), allocatable :: state_vector(:)

    state_vector = this%state_vector
  end function get_state_vector


  pure function get_physics_type(this) result(physics_type)
    class(settings_t), intent(in) :: this
    character(len=:), allocatable :: physics_type

    physics_type = trim(adjustl(this%physics_type))
  end function get_physics_type


  pure subroutine set_nb_eqs(this, nb_eqs)
    class(settings_t), intent(inout) :: this
    integer, intent(in) :: nb_eqs
    this%nb_eqs = nb_eqs
  end subroutine set_nb_eqs


  pure integer function get_nb_eqs(this)
    class(settings_t), intent(in) :: this
    get_nb_eqs = this%nb_eqs
  end function get_nb_eqs


  pure subroutine update_block_dimensions(this)
    class(settings_t), intent(inout) :: this

    call this%dims%set_block_dims(nb_eqs=this%nb_eqs, gridpts=this%grid%get_gridpts())
  end subroutine update_block_dimensions


  pure subroutine delete(this)
    class(settings_t), intent(inout) :: this

    if (allocated(this%state_vector)) deallocate(this%state_vector)
    if (allocated(this%physics_type)) deallocate(this%physics_type)
    call this%io%delete()
    call this%solvers%delete()
    call this%grid%delete()
    call this%equilibrium%delete()
  end subroutine delete

end module mod_settings
