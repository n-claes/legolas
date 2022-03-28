module mod_matrix_node
  use mod_global_variables, only: dp
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

end module mod_matrix_node