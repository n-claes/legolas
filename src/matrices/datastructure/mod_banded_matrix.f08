!> Contains types and routines to handle banded matrices. We use the same conventions
!! as explained in the LAPACK guide http://www.netlib.org/lapack/lug/node124.html.
module mod_banded_matrix
  use mod_global_variables, only: dp
  implicit none

  private

  !> type to represent a complex banded matrix
  type, public :: banded_matrix_t
    !> number of rows
    integer :: m
    !> number of columns
    integer :: n
    !> number of subdiagonals
    integer :: kl
    !> number of superdiagonals
    integer :: ku
    !> array containing the banded storage
    complex(dp), allocatable :: AB(:, :)

    contains

    procedure, public :: get_element
    procedure, public :: set_element
    procedure, public :: destroy
  end type banded_matrix_t

  public :: new_banded_matrix

contains

  !> Constructor for a new banded matrix with a given number of rows, columns,
  !! subdiagonals and superdiagonals. Allocates and initialises the datatype.
  pure function new_banded_matrix(rows, cols, subdiags, superdiags) result(matrix)
    !> number of rows
    integer, intent(in) :: rows
    !> number of columns
    integer, intent(in) :: cols
    !> number of subdiagonals
    integer, intent(in) :: subdiags
    !> number of superdiagonals
    integer, intent(in) :: superdiags
    !> banded matrix datatype
    type(banded_matrix_t) :: matrix

    matrix%m = rows
    matrix%n = cols
    matrix%kl = subdiags
    matrix%ku = superdiags
    allocate(matrix%AB(subdiags + superdiags + 1, cols))
    ! initialise all to zero
    matrix%AB = (0.0d0, 0.0d0)
  end function new_banded_matrix


  !> Retrieves the element at position (row, col) of the original matrix.
  !! See the LAPACK documentation, element $a_{ij}$ of the original matrix is stored
  !! at position $(ku + 1 + i - j, j)$ (with $ku$ the number of superdiagonals).
  pure function get_element(this, row, col) result(element)
    !> type instance
    class(banded_matrix_t), intent(in) :: this
    !> the row index of the original position
    integer, intent(in) :: row
    !> the column index of the original position
    integer, intent(in) :: col
    !> the element at original position (row, col)
    complex(dp) :: element

    if (is_within_band(matrix=this, row=row, col=col)) then
      element = this%AB(this%ku + 1 + row - col, col)
    else
      element = (0.0d0, 0.0d0)
    end if
  end function get_element


  !> Sets the element $a_{ij}$ of the original array into the banded structure.
  !! The <tt>row</tt> and <tt>col</tt> arguments refer to the row and column indices
  !! of the element in the original array. This routine has no effect if the location
  !! falls outside of the banded structure.
  pure subroutine set_element(this, row, col, element)
    !> type instance
    class(banded_matrix_t), intent(inout) :: this
    !> row index of element
    integer, intent(in) :: row
    !> column index of element
    integer, intent(in) :: col
    !> value for the element at (row, col)
    complex(dp), intent(in) :: element

    if (.not. is_within_band(matrix=this, row=row, col=col)) return
    this%AB(this%ku + 1 + row - col, col) = element
  end subroutine set_element


  !> Destructor, deallocates the datastructure.
  pure subroutine destroy(this)
    !> type instance
    class(banded_matrix_t), intent(inout) :: this

    if (allocated(this%AB)) deallocate(this%AB)
  end subroutine destroy


  !> Checks if a given position (row, col) is within the banded structure, i.e.
  !! $$ max(1, col - ku) \leq row \leq min(m, col + kl) $$
  !! with $ku$ the number of superdiagonals and $kl$ the number of subdiagonals.
  pure logical function is_within_band(matrix, row, col)
    !> the banded matrix structure
    type(banded_matrix_t), intent(in) :: matrix
    !> row index
    integer, intent(in) :: row
    !> column index
    integer, intent(in) :: col

    is_within_band = ( &
      max(1, col - matrix%ku) <= row .and. row <= min(matrix%m, col + matrix%kl) &
    )
  end function is_within_band

end module mod_banded_matrix
