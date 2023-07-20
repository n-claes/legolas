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
  contains
      procedure, public :: dspline
      procedure, public :: ddspline
      procedure, public :: delete

      procedure, private :: set_basis_function
  end type sv_component_t

contains

  function new_sv_component(name, default_spline) result(sv_comp)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: default_spline
    type(sv_component_t) :: sv_comp

    sv_comp%name = name
    call sv_comp%set_basis_function(default_spline)
  end function new_sv_component


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
