module mod_matrix_structure
  use mod_logging, only: log_message, str
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
    procedure :: delete_matrix
  end type matrix_T

  public :: new_matrix

contains

  pure function new_matrix(nb_rows) result(matrix)
    integer, intent(in) :: nb_rows
    type(matrix_t) :: matrix
    integer :: i


    matrix%matrix_dim = nb_rows
    allocate(matrix%rows(nb_rows))
    do i = 1, matrix%matrix_dim
      matrix%rows(i) = new_row()
    end do
  end function new_matrix


  subroutine add_element(this, row, column, element)
    class(matrix_t), intent(inout) :: this
    integer, intent(in) :: row
    integer, intent(in) :: column
    class(*), intent(in) :: element

    if (.not. is_valid_element(element)) return
    if (.not. is_valid_index(this, row)) return
    if (.not. is_valid_index(this, column)) return

    call this%rows(row)%add_node(column, element)
  end subroutine add_element


  logical function is_valid_element(element) result(is_valid)
    use mod_global_variables, only: dp
    use mod_check_values, only: is_zero

    class(*), intent(in) :: element

    is_valid = .false.
    select type(element)
      type is (real(dp))
        is_valid = (.not. is_zero(element))
      type is (complex(dp))
        is_valid = (.not. is_zero(element))
      class default
        call log_message("is_valid_element: unexpected element type", level="error")
    end select
  end function is_valid_element


  logical function is_valid_index(matrix, index) result(is_valid)
    type(matrix_t), intent(in) :: matrix
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


  pure subroutine delete_matrix(this)
    class(matrix_t), intent(inout) :: this
    integer :: i

    do i = 1, this%matrix_dim
      call this%rows(i)%delete_row()
    end do
    deallocate(this%rows)
  end subroutine delete_matrix

end module mod_matrix_structure