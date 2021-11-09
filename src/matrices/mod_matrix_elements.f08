!=================================================================
!> Module defining a linked list implementation to handle addition of matrix elements.
module mod_matrix_elements
  use mod_global_variables, only: dp, str_len
  use mod_logging, only: log_message, str
  implicit none

  private

  !> node in the linked list
  type :: node_t
    !> the matrix element
    complex(dp) :: matrix_element
    !> the state vector index position of the corresponding element
    integer :: position(2)
    !> pointer to the next node
    type(node_t), pointer :: next
  end type


  !> linked list containing the matrix elements and corresponding positions
  type, public :: matrix_elements_t
    !> pointer to the current head, which is the last added element
    type(node_t), pointer, private :: head
    !> if the linked list is initialised or not
    logical, private :: initialised = .false.
    !> the total number of elements added
    integer, private :: nb_elements

    contains

    procedure, private :: add_real
    procedure, private :: add_complex
    procedure, private :: initialise

    generic :: add => add_real, add_complex
    procedure :: get_values
    procedure :: get_positions
    procedure :: get_size
    procedure :: delete
  end type

contains


  !> Adds a real element with a given location to the list. The element is converted
  !! to a complex quantity and `add_complex` is called.
  subroutine add_real(this, element, location)
    !> type instance
    class(matrix_elements_t) :: this
    !> the element to add
    real(dp), intent(in) :: element
    !> name of state vector locations
    character(len=*), intent(in) :: location(2)

    call this%add_complex(element=cmplx(element, kind=dp), location=location)
  end subroutine add_real


  !> Adds a complex element with a given location to the list. If an element is added
  !! with a location that does not correspond to a state vector variable name
  !! the element is skipped and not added to the list.
  subroutine add_complex(this, element, location)
    use mod_global_variables, only: state_vector
    use mod_get_indices, only: get_index

    !> type instance
    class(matrix_elements_t) :: this
    !> the element to add
    complex(dp), intent(in) :: element
    !> name of state vector locations
    character(len=*), intent(in) :: location(2)

    type(node_t), pointer :: node
    integer :: position(2)

    position = get_index(location, state_vector)
    if (.not. is_valid_position(position)) return

    if (.not. this%initialised) call this%initialise()

    allocate(node)
    node%matrix_element = element
    node%position = position
    node%next => this%head
    this%head => node
    this%nb_elements = this%nb_elements + 1
  end subroutine add_complex


  !> Initialises the linked list
  subroutine initialise(this)
    !> type instance
    class(matrix_elements_t) :: this

    this%head => null()
    this%initialised = .true.
    this%nb_elements = 0
  end subroutine initialise


  !> Checks if a given position is valid. Only returns `.true.` if both position
  !! indices are non-zero.
  logical function is_valid_position(position)
    integer, intent(in) :: position(:)

    is_valid_position = (all(position > 0))
  end function is_valid_position


  !> Returns the matrix elements in the list.
  function get_values(this) result(values)
    !> type instance
    class(matrix_elements_t) :: this
    !> matrix elements in the list
    complex(dp) :: values(this%nb_elements)
    !> pointer to the current head
    type(node_t), pointer :: current_head
    integer :: i

    current_head => this%head
    do i = this%nb_elements, 1, -1
      values(i) = current_head%matrix_element
      current_head => current_head%next
    end do
  end function get_values


  !> Returns the position indices of the matrix elements in the list.
  function get_positions(this) result(positions)
    !> type instance
    class(matrix_elements_t) :: this
    !> position indices of matrix elements in the list
    integer :: positions(this%nb_elements, 2)
    !> pointer to the current head
    type(node_t), pointer :: current_head
    integer :: i

    current_head => this%head
    do i = this%nb_elements, 1, -1
      positions(i, :) = current_head%position
      current_head => current_head%next
    end do
  end function get_positions


  !> Returns the number of elements in the list.
  integer function get_size(this)
    !> type instance
    class(matrix_elements_t) :: this

    get_size = this%nb_elements
  end function get_size


  !> Deletes the list with matrix elements, deallocates the pointers in the
  !! linked list by recursively descending starting from the top head.
  subroutine delete(this)
    !> type instance
    class(matrix_elements_t) :: this

    call delete_node(this%head)
    this%initialised = .false.

    contains

    !> Deletes a single node in the linked list
    recursive subroutine delete_node(head)
      !> the current head that gets deallocated
      type(node_t), pointer, intent(inout) :: head
      !> pointer to the current head to pass recursively
      type(node_t), pointer :: current_head

      current_head => head
      deallocate(head)
      if (associated(current_head%next)) call delete_node(current_head%next)
    end subroutine delete_node
  end subroutine delete

end module mod_matrix_elements