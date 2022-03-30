module mod_matrix_structure
  use mod_global_variables, only: dp
  use mod_logging, only: log_message, str
  use mod_matrix_node, only: node_t
  use mod_matrix_row, only: row_t, new_row
  implicit none

  private

  !> General matrix type, represents the linked list implementation
  type, public :: matrix_t
    !> dimension of the matrix, number of rows
    integer :: matrix_dim
    !> array containing the various rows
    type(row_t), allocatable :: rows(:)

    contains

    procedure :: add_element
    procedure :: get_real_element
    procedure :: get_complex_element
    procedure :: get_total_nb_elements
    procedure :: delete_matrix

    procedure, private :: add_matrices
    generic :: operator(+) => add_matrices
  end type matrix_T

  ! interface for submodule implementations
  interface
    !> Overloads the addition operator between two matrix datastructures
    module function add_matrices(matrix1, matrix2) result(matrix)
      !> left-hand side of + operator
      class(matrix_t), intent(in) :: matrix1
      !> right-hand side of + operator
      class(matrix_t), intent(in) :: matrix2
      !> result of addition
      type(matrix_t) :: matrix
    end function add_matrices
  end interface

  public :: new_matrix

contains

  !> Constructor for a new matrix matrix with a given number of rows.
  !! Allocates and initialises the matrix datatype.
  pure function new_matrix(nb_rows) result(matrix)
    !> number of rows in the matrix
    integer, intent(in) :: nb_rows
    !> matrix datatype with rows/columns in a linked list
    type(matrix_t) :: matrix
    integer :: i

    matrix%matrix_dim = nb_rows
    allocate(matrix%rows(nb_rows))
    do i = 1, matrix%matrix_dim
      matrix%rows(i) = new_row()
    end do
  end function new_matrix


  !> Adds a given element at a certain (row, column) position to the matrix
  !! datastructure. Elements that are zero are not added, sanity checks are done
  !! on the row and column indices.
  subroutine add_element(this, row, column, element)
    !> type instance
    class(matrix_t), intent(inout) :: this
    !> row position of the element
    integer, intent(in) :: row
    !> column position of the element
    integer, intent(in) :: column
    !> polymorphic variable to add to the matrix
    class(*), intent(in) :: element

    if (.not. is_valid_element(element)) return
    if (.not. is_valid_index(this, row)) return
    if (.not. is_valid_index(this, column)) return

    call this%rows(row)%add_node(column, element)
  end subroutine add_element


  !> Checks if a given element is valid in order to add it to the matrix.
  !! Returns `.true.` if the element is of type real or complex, `.false.` otherwise.
  logical function is_valid_element(element) result(is_valid)
    use mod_check_values, only: is_zero

    !> Matrix element that is to be added
    class(*), intent(in) :: element

    is_valid = .false.
    select type(element)
      type is (real(dp))
        is_valid = (.not. is_zero(element))
      type is (complex(dp))
        is_valid = (.not. is_zero(element))
      class default
        call log_message("adding unexpected element type", level="error")
    end select
  end function is_valid_element


  !> Checks if a given index is valid for the current matrix datastructure.
  !! Returns `.true.` if the index (either row or column) is larger than 0 and
  !! smaller than the dimension of the matrix. Returns `.false.` otherwise.
  logical function is_valid_index(matrix, index) result(is_valid)
    !> matrix datastructure object
    type(matrix_t), intent(in) :: matrix
    !> index to check
    integer, intent(in) :: index

    is_valid = .true.
    if (index <= 0 .or. index > matrix%matrix_dim) then
      call log_message( &
        "row/column index " // str(index) // " is outside of matrix dimension", &
        level="error" &
      )
      is_valid = .false.
    end if
  end function is_valid_index


  !> Returns the real element associated with the linked-list node at position
  !! (row, column) in the matrix datastructure. Throws appropriate errors if the
  !! node element types do not match. Non-existing nodes correspond to zero elements,
  !! so in that case this function returns (real) zero.
  function get_real_element(this, row, column) result(element)
    !> type instance
    class(matrix_t), intent(in) :: this
    !> row position of the needed element
    integer, intent(in) :: row
    !> column position of the needed element
    integer, intent(in) :: column
    !> the element at position (row, column) in the matrix
    real(dp) :: element
    type(node_t), pointer :: node

    ! if node is not present then element is zero
    element = 0.0d0
    node => this%rows(row)%get_node(column=column)
    if (associated(node)) call node%get_node_element(element)
    nullify(node)
  end function get_real_element


  !> Returns the complex element associated with the linked-list node at position
  !! (row, column) in the matrix datastructure. Throws appropriate errors if the
  !! node element types do not match. Non-existing nodes correspond to zero elements,
  !! so in that case this function returns (complex) zero.
  function get_complex_element(this, row, column) result(element)
    !> type instance
    class(matrix_t), intent(in) :: this
    !> row position of the needed element
    integer, intent(in) :: row
    !> column position of the needed element
    integer, intent(in) :: column
    !> the element at position (row, column) in the matrix
    complex(dp) :: element
    type(node_t), pointer :: node

    element = (0.0d0, 0.0d0)
    node => this%rows(row)%get_node(column=column)
    if (associated(node)) call node%get_node_element(element)
    nullify(node)
  end function get_complex_element


  !> Returns the total number of elements (nodes) across the various rows.
  pure function get_total_nb_elements(this) result(total_nb_elements)
    !> type instance
    class(matrix_t), intent(in) :: this
    !> total number of (non-zero) elements in this matrix
    integer :: total_nb_elements
    integer :: i

    total_nb_elements = 0
    do i = 1, this%matrix_dim
      total_nb_elements = total_nb_elements + this%rows(i)%nb_elements
    end do
  end function get_total_nb_elements


  !> Deallocates the matrix datastructure, nullifies all corresponding pointers and
  !! deallocates the various nodes in the rows.
  pure subroutine delete_matrix(this)
    class(matrix_t), intent(inout) :: this
    integer :: i

    if (.not. allocated(this%rows)) return
    do i = 1, this%matrix_dim
      call this%rows(i)%delete_row()
    end do
    deallocate(this%rows)
  end subroutine delete_matrix

end module mod_matrix_structure