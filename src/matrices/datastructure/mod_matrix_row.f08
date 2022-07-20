! =============================================================================
!> Module that contains the implementation of a single row of nodes in the
!! linked-list matrix representation.
module mod_matrix_row
  use mod_logging, only: log_message, str
  use mod_matrix_node, only: node_t, new_node
  implicit none

  private

  !> Linked list for a given row index, contains the column values
  type, public :: row_t
    !> number of elements in linked list
    integer :: nb_elements
    !> pointer to head (first element added)
    type(node_t), pointer :: head
    !> pointer to tail (last element added)
    type(node_t), pointer :: tail

    contains

    procedure :: add_node
    procedure :: get_node
    procedure :: delete_row
    procedure :: delete_node_from_row

    procedure, private :: create_first_node
    procedure, private :: append_node
  end type row_t

  public :: new_row

contains

  !> Constructor for a new row, initialises the linked list datastructure
  !! and sets the current head and tail pointers to `null()`.
  pure function new_row() result(row)
    type(row_t) :: row

    row%nb_elements = 0
    row%head => null()
    row%tail => null()
  end function new_row


  !> Adds a new node to the linked list with a given column index and value.
  subroutine add_node(this, column, element)
    !> type instance
    class(row_t), intent(inout) :: this
    !> column position of the element
    integer, intent(in) :: column
    !> the element to be added
    class(*), intent(in) :: element

    if (.not. associated(this%head)) then
      call this%create_first_node(column, element)
    else
      call this%append_node(column, element)
    end if
  end subroutine add_node


  !> Subroutine to add the first node to the linked list. Allocates a new node and
  !! sets both the head and tail to this node.
  pure subroutine create_first_node(this, column, element)
    !> type instance
    class(row_t), intent(inout) :: this
    !> column position of element
    integer, intent(in) :: column
    !> the element to be added
    class(*), intent(in) :: element

    allocate(this%head, source=new_node(column, element))
    this%tail => this%head
    this%nb_elements = this%nb_elements + 1
  end subroutine create_first_node


  !> Subroutine to append a new node to an already existing list of nodes.
  !! A new node is created, appended, and the tail is updated.
  subroutine append_node(this, column, element)
    !> type instance
    class(row_t), intent(inout) :: this
    !> column position of element
    integer, intent(in) :: column
    !> the element to be added
    class(*), intent(in) :: element

    type(node_t), pointer :: node

    node => this%get_node(column)
    ! check if node already exists
    if (associated(node)) then
      call node%add_to_node_element(element)
    else
      allocate(this%tail%next, source=new_node(column, element))
      ! update tail to last element added
      this%tail => this%tail%next
      this%nb_elements = this%nb_elements + 1
    end if
    nullify(node)
  end subroutine append_node


  !> Returns a pointer to the node corresponding to the given column.
  !! Returns a nullified pointer if no node containing the given column index
  !! was found.
  function get_node(this, column) result(node)
    !> type instance
    class(row_t), intent(in) :: this
    !> column index
    integer, intent(in) :: column
    !> the node with a column value that matches column
    type(node_t), pointer :: node

    type(node_t), pointer :: current_node
    integer :: i

    node => null()
    current_node => this%head
    ! loop over nodes, return node if column index matches
    do i = 1, this%nb_elements
      if (column == current_node%column) then
        node => current_node
        nullify(current_node)
        exit
      end if
      current_node => current_node%next
    end do
    nullify(current_node)
  end function get_node


  !> Deletes a given node from the current row.
  subroutine delete_node_from_row(this, column)
    !> type instance
    class(row_t), intent(inout) :: this
    !> column index of node to be deleted
    integer, intent(in) :: column

    type(node_t), pointer :: node
    type(node_t), pointer :: next_node
    integer :: i

    ! do nothing for empty rows
    if (.not. associated(this%head)) return

    node => this%head
    ! check if head is the one being deleted
    if (column == node%column) then
      this%head => this%head%next
      call node%delete()
      call decrement_and_nullify()
      return
    end if

    next_node => this%head%next
    ! -1 since we're working with current and next nodes, head and tail are checked
    do i = 1, this%nb_elements - 1
      ! check if next node is tail and will be deleted
      if (associated(next_node, this%tail) .and. column == next_node%column) then
        call next_node%delete()
        node%next => null()
        this%tail => node
        call decrement_and_nullify()
        exit
      end if
      if (column == next_node%column) then
        node%next => next_node%next
        call next_node%delete()
        call decrement_and_nullify()
        exit
      end if
      node => next_node
      next_node => next_node%next
    end do
    nullify(node)
    nullify(next_node)

    contains

    !> Decrements the total number of nodes in the row and nullifies the node pointers
    !! of the containing subroutine. Called whenever a node gets deleted.
    subroutine decrement_and_nullify()
      this%nb_elements = this%nb_elements - 1
      nullify(node)
      nullify(next_node)
    end subroutine decrement_and_nullify
  end subroutine delete_node_from_row


  !> Deletes a given linked list row by recursively iterating over all nodes.
  !! Nullifies the pointers and deallocates the elements.
  pure subroutine delete_row(this)
    !> type instance
    class(row_t), intent(inout) :: this

    if (associated(this%head)) call delete_node(node=this%head)
    nullify(this%head)
    nullify(this%tail)
    this%nb_elements = 0

    contains

    !> Recursive subroutine that deallocates a given node in the linked list.
    pure recursive subroutine delete_node(node)
      type(node_t), intent(inout) :: node
      type(node_t), pointer :: next_node

      next_node => null()
      if (associated(node%next)) next_node => node%next
      call node%delete()
      if (associated(next_node)) call delete_node(next_node)
      nullify(next_node)
    end subroutine delete_node
  end subroutine delete_row
end module mod_matrix_row
