module mod_matrix_node
  use mod_global_variables, only: dp, NaN
  use mod_logging, only: log_message
  implicit none

  private

  !> Base node corresponding to a given column - value pair for a given row
  type, public :: node_t
    !> column index
    integer :: column
    !> value for the matrix element, can be real/complex
    class(*), allocatable :: element
    !> pointer to next node
    class(node_t), pointer :: next

    contains

    procedure, private :: get_real_node_element
    generic :: get_node_element => get_real_node_element
  end type node_t

  public :: new_node

contains

  pure function new_node(column, element) result(node)
    integer, intent(in) :: column
    class(*), intent(in) :: element
    type(node_t) :: node

    allocate(node%element, source=element)
    node%column = column
    node%next => null()
  end function new_node


  subroutine get_real_node_element(this, element)
    class(node_t), intent(in) :: this
    real(dp), intent(out) :: element

    select type(item => this%element)
      type is (real(dp))
        element = item
      class default
        element = NaN
        call throw_type_error(element_type="real")
    end select

  end subroutine get_real_node_element


  subroutine throw_type_error(element_type)
    character(len=*), intent(in) :: element_type

    call log_message("node element does not have type " // element_type, level="error")
  end subroutine throw_type_error

end module mod_matrix_node