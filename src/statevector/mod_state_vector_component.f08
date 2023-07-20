module mod_state_vector_component
  use mod_global_variables, only: str_len_arr
  use mod_basis_functions, only: basis_function, hquad, dhquad, hcubic, dhcubic, &
    ddhcubic
  use mod_logging, only: logger
  implicit none

  private

  type, public :: sv_component_t
    procedure(basis_function), pointer, nopass :: spline => null()
    character(len=str_len_arr), private :: name
    logical, private :: basis_function_is_set = .false.
  contains
      procedure, public :: set_basis_function
      procedure, public :: has_basis_function
      procedure, public :: delete
  end type sv_component_t

contains

  function new_sv_component(name, default_spline) result(sv_comp)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: default_spline
    type(sv_component_t) :: sv_comp

    sv_comp%name = name
    sv_comp%basis_function_is_set = .false.
    call sv_comp%set_basis_function(default_spline)
  end function new_sv_component


  subroutine set_basis_function(this, spline_name)
    use mod_basis_function_names

    class(sv_component_t), intent(inout) :: this
    character(len=*), intent(in) :: spline_name

    select case(spline_name)
    case(QUADRATIC)
      this%spline => hquad
    case(DQUADRATIC)
      this%spline => dhquad
    case(CUBIC)
      this%spline => hcubic
    case(DCUBIC)
      this%spline => dhcubic
    case(DDCUBIC)
      this%spline => ddhcubic
    case default
      call logger%error("unknown basis function name: " // trim(adjustl(spline_name)))
      return
    end select
    this%basis_function_is_set = .true.
  end subroutine set_basis_function


  pure logical function has_basis_function(this)
    class(sv_component_t), intent(in) :: this
    has_basis_function = this%basis_function_is_set
  end function has_basis_function


  pure subroutine delete(this)
    class(sv_component_t), intent(inout) :: this
    this%spline => null()
  end subroutine delete

end module mod_state_vector_component
