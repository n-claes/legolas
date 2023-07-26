module mod_state_vector_component
  use mod_global_variables, only: str_len_arr
  use mod_basis_function_names, only: QUADRATIC, CUBIC
  use mod_basis_functions, only: basis_function
  use mod_logging, only: logger
  implicit none

  private

  type, public :: sv_component_t
    procedure(basis_function), pointer, nopass :: spline => null()
    character(len=str_len_arr), private :: spline_name
    character(len=str_len_arr), private :: name
    logical, private :: is_initialised = .false.
  contains
    procedure, public :: dspline
    procedure, public :: ddspline

    procedure, public :: get_name
    procedure, public :: get_basis_function_name
    procedure, public :: get_default_basis_function
    procedure, public :: set_basis_function
    procedure, public :: delete
  end type sv_component_t

  public :: new_sv_component

contains

  function new_sv_component(name) result(sv_comp)
    character(len=*), intent(in) :: name
    type(sv_component_t) :: sv_comp

    sv_comp%name = name
    sv_comp%spline_name = ""
    sv_comp%is_initialised = .true.
  end function new_sv_component


  function get_name(this) result(name)
    class(sv_component_t), intent(in) :: this
    character(:), allocatable :: name
    name = this%name
  end function get_name


  function get_basis_function_name(this) result(name)
    class(sv_component_t), intent(in) :: this
    character(:), allocatable :: name
    name = this%spline_name
  end function get_basis_function_name


  function get_default_basis_function(this) result(name)
    use mod_state_vector_names

    class(sv_component_t), intent(in) :: this
    character(:), allocatable :: name

    select case(this%name)
    case(sv_rho1_name, sv_v2_name, sv_v3_name, sv_T1_name, sv_a1_name)
      name = QUADRATIC
    case(sv_v1_name, sv_a2_name, sv_a3_name)
      name = CUBIC
    case default
      call logger%error( &
        "default basis function not defined for state vector component " &
        // trim(adjustl(this%name)) &
      )
    end select
  end function get_default_basis_function


  subroutine set_basis_function(this, spline_name)
    use mod_basis_functions, only: hquad, hcubic

    class(sv_component_t), intent(inout) :: this
    character(len=*), intent(in) :: spline_name

    this%spline_name = spline_name
    select case(this%spline_name)
      case(QUADRATIC)
        this%spline => hquad
      case(CUBIC)
        this%spline => hcubic
      case default
        call logger%error("unknown basis function name: " // trim(adjustl(spline_name)))
        return
    end select
  end subroutine set_basis_function


  function dspline(this) result(spline)
    use mod_basis_functions, only: dhquad, dhcubic

    class(sv_component_t), intent(in) :: this
    procedure(basis_function), pointer :: spline

    spline => null()
    select case(this%spline_name)
      case(QUADRATIC)
        spline => dhquad
      case(CUBIC)
        spline => dhcubic
      case default
        call logger%error( &
          trim(adjustl(this%spline_name)) &
          // " basis function has no defined derivative" &
        )
        return
    end select
  end function dspline


  function ddspline(this) result(spline)
    use mod_basis_functions, only: ddhcubic

    class(sv_component_t), intent(in) :: this
    procedure(basis_function), pointer :: spline

    spline => null()
    select case(this%spline_name)
      case(CUBIC)
        spline => ddhcubic
      case default
        call logger%error( &
          trim(adjustl(this%spline_name)) &
          // " basis function has no defined second derivative" &
        )
        return
    end select
  end function ddspline


  pure subroutine delete(this)
    class(sv_component_t), intent(inout) :: this
    this%spline => null()
  end subroutine delete

end module mod_state_vector_component
