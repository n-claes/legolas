module mod_matrix_row
  use mod_logging, only: log_message, str
  use mod_matrix_node, only: node_t
  implicit none

  private

  !> Linked list for a given row index, contains the column values
  type, public :: row_t
    !> initialisation check
    logical :: initialised = .false.
    !> number of elements in linked list
    integer :: nb_elements
    !> pointer to current head (last added element)
    type(node_t), pointer :: head

    contains

    procedure :: add_node
    procedure :: get_node
    procedure :: delete_row
  end type row_t

  public :: new_row

contains

  !> Constructor for a new row, initialises the linked list datastructure
  !! and sets the current head pointer to `null()`.
  pure function new_row() result(row)
    type(row_t) :: row

    row%initialised = .true.
    row%nb_elements = 0
    row%head => null()
  end function new_row


  !> Adds a new node to the linked list with a given column index and value.
  pure subroutine add_node(this, column, element)
    use mod_matrix_node, only: new_node

    !> type instance
    class(row_t), intent(inout) :: this
    !> column position of the element
    integer, intent(in) :: column
    !> the element to be added
    class(*), intent(in) :: element

    if (.not. associated(this%head)) then
      ! list is empty, allocate new node and set as head
      allocate(this%head, source=new_node(column, element))
    else
      ! list is not empty, allocate new node and set as next
      allocate(this%head%next, source=new_node(column, element))
      ! update head to last added node
      this%head => this%head%next
    end if
    this%nb_elements = this%nb_elements + 1
  end subroutine add_node


  !> Returns the node corresponding to the given column.
  !! Throws an error if no node containing the given column index was found.
  function get_node(this, column) result(node)
    !> type instance
    class(row_t), intent(in) :: this
    !> column index
    integer, intent(in) :: column
    !> the node with a column value that matches column
    type(node_t) :: node
    type(node_t), pointer :: current_node
    integer :: i

    current_node => this%head
    do i = 1, this%nb_elements
      if (column == current_node%column) then
        node = current_node
        current_node => null()
        exit
      end if
      current_node => this%head%next
    end do
    current_node => null()
    ! if node has not been found then element is still deallocated
    if (.not. allocated(node%element)) then
      call log_message( &
        "node with column index " // str(column) // " was not found", &
        level="error" &
      )
    end if
  end function get_node


  !> Deletes a given linked list row by recursively iterating over all nodes.
  !! Nullifies the pointers and deallocates the elements.
  pure subroutine delete_row(this)
    class(row_t), intent(inout) :: this

    if (associated(this%head)) call delete_node(node=this%head)
    this%initialised = .false.
    this%nb_elements = 0

    contains

    !> Recursive subroutine that deallocates a given node in the linked list.
    pure recursive subroutine delete_node(node)
      type(node_t), intent(inout) :: node
      type(node_t), pointer :: next_node

      next_node => null()
      if (associated(node%next)) next_node => node%next
      deallocate(node%element)
      node%next => null()
      if (associated(next_node)) call delete_node(next_node)
      next_node => null()
    end subroutine delete_node
  end subroutine delete_row
end module mod_matrix_row