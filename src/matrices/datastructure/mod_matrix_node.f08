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

    procedure :: add_to_node_element
    procedure, private :: add_real_to_node_element

    procedure, private :: get_real_node_element
    generic :: get_node_element => get_real_node_element
  end type node_t

  public :: new_node

contains

  !> Constructor for a new node, sets the column and element attributes.
  !! No nodes are linked yet; the pointer to the next node is initialised to `null()`.
  pure function new_node(column, element) result(node)
    !> column index
    integer, intent(in) :: column
    !> element added to the node
    class(*), intent(in) :: element
    !> new node with given column and element attributes
    type(node_t) :: node

    allocate(node%element, source=element)
    node%column = column
    node%next => null()
  end function new_node


  !> Generic type-checking function to add a polymorphic element to the existing
  !! element attribute.
  subroutine add_to_node_element(this, element)
    !> type instance
    class(node_t), intent(inout) :: this
    !> element to add to the existing element
    class(*), intent(in) :: element

    select type(element)
      type is (real(dp))
        call this%add_real_to_node_element(element)
    end select
  end subroutine add_to_node_element


  !> Sums the current real node element with a new element. The "old" element is
  !! deallocated and reallocated with its new value.
  subroutine add_real_to_node_element(this, element_to_add)
    !> type instance
    class(node_t), intent(inout) :: this
    !> element to add to the existing element
    real(dp), intent(in) :: element_to_add
    real(dp) :: existing_element

    call this%get_node_element(existing_element)
    deallocate(this%element)
    allocate(this%element, source=(existing_element + element_to_add))
  end subroutine add_real_to_node_element


  !> Getter for nodes with real elements, returns a real `element` attribute.
  !! Throws an error if the element attribute is not of type real.
  subroutine get_real_node_element(this, element)
    !> type instance
    class(node_t), intent(in) :: this
    !> corresponding element
    real(dp), intent(out) :: element

    select type(item => this%element)
      type is (real(dp))
        element = item
      class default
        element = NaN
        call throw_type_error(element_type="real")
    end select
  end subroutine get_real_node_element


  !> Throws an error message stating that the requested element type and the
  !! corresponding element type do not match.
  subroutine throw_type_error(element_type)
    !> element type to display in the error message
    character(len=*), intent(in) :: element_type

    call log_message("node element does not have type " // element_type, level="error")
  end subroutine throw_type_error

end module mod_matrix_node