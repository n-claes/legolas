module mod_state_vector
  use mod_state_vector_component, only: sv_component_t
  use mod_logging, only: logger, str
  implicit none

  private

  type, public :: state_vector_t
    type(sv_component_t), allocatable :: components(:)

  contains
    procedure, public :: assemble
    procedure, public :: set_basis_functions
    procedure, public :: get_names
    procedure, public :: get_basis_functions
    procedure, public :: delete

    procedure, private :: set_default_basis_functions
    procedure, private :: is_compatible_with
  end type state_vector_t

  logical, save, private :: sv_components_initialised = .false.
  type(sv_component_t), public, protected :: sv_rho1
  type(sv_component_t), public, protected :: sv_v1
  type(sv_component_t), public, protected :: sv_v2
  type(sv_component_t), public, protected :: sv_v3
  type(sv_component_t), public, protected :: sv_T1
  type(sv_component_t), public, protected :: sv_a1
  type(sv_component_t), public, protected :: sv_a2
  type(sv_component_t), public, protected :: sv_a3

contains

  subroutine assemble(this, physics_type)
    class(state_vector_t), intent(inout) :: this
    character(len=*), intent(in) :: physics_type

    if (.not. sv_components_initialised) call initialise_sv_components()
    select case(physics_type)
      case("hd")
        this%components = [sv_rho1, sv_v1, sv_v2, sv_v3, sv_T1]
      case("hd-1d")
        this%components = [sv_rho1, sv_v1, sv_T1]
      case default
        this%components = [sv_rho1, sv_v1, sv_v2, sv_v3, sv_T1, sv_a1, sv_a2, sv_a3]
    end select
  end subroutine assemble


  subroutine set_basis_functions(this, splines)
    class(state_vector_t), intent(inout) :: this
    character(len=*), intent(in) :: splines(:)
    integer :: i

    if (size(splines) == 0) then
      call this%set_default_basis_functions()
      return
    end if
    if (.not. this%is_compatible_with(splines)) return
    do i = 1, size(this%components)
      call this%components(i)%set_basis_function(splines(i))
    end do
  end subroutine set_basis_functions


  subroutine set_default_basis_functions(this)
    class(state_vector_t), intent(inout) :: this
    character(:), allocatable :: default_name
    integer :: i

    do i = 1, size(this%components)
      default_name = this%components(i)%get_default_basis_function()
      call this%components(i)%set_basis_function(default_name)
    end do
  end subroutine set_default_basis_functions


  function get_names(this) result(names)
    class(state_vector_t), intent(in) :: this
    character(len=:), allocatable :: names(:)
    integer :: i

    allocate(names(size(this%components)), mold=this%components(1)%get_name())
    do i = 1, size(this%components)
      names(i) = this%components(i)%get_name()
    end do
  end function get_names


  function get_basis_functions(this) result(names)
    class(state_vector_t), intent(in) :: this
    character(len=:), allocatable :: names(:)
    integer :: i

    allocate(names(size(this%components)), mold=this%components(1)%get_name())
    do i = 1, size(this%components)
      names(i) = this%components(i)%get_basis_function_name()
    end do
  end function get_basis_functions


  logical function is_compatible_with(this, splines)
    class(state_vector_t), intent(in) :: this
    character(len=*), intent(in) :: splines(:)
    integer :: size_sv, size_splines

    size_sv = size(this%components)
    size_splines = size(splines)
    is_compatible_with = .true.
    if (size_sv == size_splines) return

    is_compatible_with = .false.
    call logger%error( &
      "state vector size (" // str(size_sv) // ") is not compatible with " &
      // "given number of basis functions (" // str(size_splines) // ")" &
    )
  end function is_compatible_with


  subroutine delete(this)
    class(state_vector_t), intent(inout) :: this

    if (allocated(this%components)) deallocate(this%components)
    call sv_rho1%delete()
    call sv_v1%delete()
    call sv_v2%delete()
    call sv_v3%delete()
    call sv_T1%delete()
    call sv_a1%delete()
    call sv_a2%delete()
    call sv_a3%delete()
    sv_components_initialised = .false.
  end subroutine delete


  subroutine initialise_sv_components()
    use mod_state_vector_component, only: new_sv_component
    use mod_state_vector_names

    sv_rho1 = new_sv_component(sv_rho1_name)
    sv_v1 = new_sv_component(sv_v1_name)
    sv_v2 = new_sv_component(sv_v2_name)
    sv_v3 = new_sv_component(sv_v3_name)
    sv_T1 = new_sv_component(sv_T1_name)
    sv_a1 = new_sv_component(sv_a1_name)
    sv_a2 = new_sv_component(sv_a2_name)
    sv_a3 = new_sv_component(sv_a3_name)
    sv_components_initialised = .true.
  end subroutine initialise_sv_components


end module mod_state_vector
