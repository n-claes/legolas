!> Contains types and routines to handle banded Hermitian matrices.
!! We use the same conventions as explained in the LAPACK guide
!! http://www.netlib.org/lapack/lug/node124.html.
module mod_banded_matrix_hermitian
  use mod_global_variables, only: dp
  use mod_logging, only: log_message, str
  implicit none

  private

  !> type to represent a complex Hermitian banded matrix
  type, public :: hermitian_banded_matrix_t
    !> number of rows/columns
    integer :: n
    !> store upper ("U") or lower ("L") triangular part
    character :: uplo
    !> number of sub/superdiagonals
    integer :: kd
    !> the matrix in banded storage
    complex(dp), allocatable :: AB(:, :)

    contains

    procedure, public :: get_element
    procedure, public :: set_element
    procedure, public :: get_total_nb_elements
    procedure, public :: is_compatible_with
    procedure, public :: destroy
  end type hermitian_banded_matrix_t

  public :: new_hermitian_banded_matrix

contains

  !> Constructor for a new Hermitian banded matrix with a given number of rows and
  !! diagonals. Allocates and initialises the datatype.
  function new_hermitian_banded_matrix(rows, diags, uplo) result(matrix)
    !> number of rows/columns
    integer, intent(in) :: rows
    !> number of sub/superdiagonals
    integer, intent(in) :: diags
    !> store upper ("U") or lower ("L") triangular part
    character, intent(in) :: uplo
    !> the Hermitian matrix in banded storage
    type(hermitian_banded_matrix_t) :: matrix

    if (.not. uplo_is_valid(uplo)) then
      call log_message( &
        "invalid uplo argument, expected 'U' or 'L', got '" // uplo // "'", &
        level="error" &
      )
      return
    end if

    matrix%n = rows
    matrix%kd = diags
    matrix%uplo = uplo
    allocate(matrix%AB(diags + 1, rows))
    !> initialise all to zero
    matrix%AB = (0.0_dp, 0.0_dp)
  end function new_hermitian_banded_matrix


  !> Retrieves the element at position (row, col) of the original matrix.
  !! See the LAPACK documentation, element $a_{ij}$ of the original matrix is stored
  !! at position $(kd + 1 + i - j, j)$ if `uplo = "U"` and at position
  !! $(1 + i - j, j)$ if `uplo = "L"` in the banded storage.
  pure function get_element(this, row, col) result(element)
    !> type instance
    class(hermitian_banded_matrix_t), intent(in) :: this
    !> the row index of the original position
    integer, intent(in) :: row
    !> the column index of the original position
    integer, intent(in) :: col
    !> the element at the original position (row, col)
    complex(dp) :: element

    if (this%uplo == "U") then
      if (is_within_band(matrix=this, row=row, col=col)) then
        element = this%AB(this%kd + 1 + row - col, col)
        ! check if transpose element is within band, if so return transpose conjugate
      else if (is_within_band(matrix=this, row=col, col=row)) then
        element = conjg(this%AB(this%kd + 1 + col - row, row))
      else
        element = (0.0_dp, 0.0_dp)
      end if
    else  ! uplo == "L"
      if (is_within_band(matrix=this, row=row, col=col)) then
        element = this%AB(1 + row - col, col)
      else if (is_within_band(matrix=this, row=col, col=row)) then
        element = conjg(this%AB(1 + col - row, row))
      else
        element = (0.0_dp, 0.0_dp)
      end if
    end if
  end function get_element


  !> Sets the element $a_{ij}$ of the original array into the banded structure.
  !! The <tt>row</tt> and <tt>col</tt> arguments refer to the row and column indices
  !! of the element in the original array. This routine has no effect if the location
  !! falls outside of the banded structure.
  pure subroutine set_element(this, row, col, element)
    !> type instance
    class(hermitian_banded_matrix_t), intent(inout) :: this
    !> row index of element
    integer, intent(in) :: row
    !> column index of element
    integer, intent(in) :: col
    !> value for the element at (row, col)
    complex(dp), intent(in) :: element

    if (.not. is_within_band(matrix=this, row=row, col=col)) return
    if (this%uplo == "U") then
      this%AB(this%kd + 1 + row - col, col) = element
    else  ! this%uplo == "L"
      this%AB(1 + row - col, col) = element
    end if
  end subroutine set_element


  !> Destructor, deallocates the datastructure.
  pure subroutine destroy(this)
    !> type instance
    class(hermitian_banded_matrix_t), intent(inout) :: this

    if (allocated(this%AB)) deallocate(this%AB)
  end subroutine destroy


  !> Returns the total number of elements inside the banded matrix
  pure integer function get_total_nb_elements(this)
    !> type instance
    class(hermitian_banded_matrix_t), intent(in) :: this
    integer :: i

    ! N elements on diagonal, square matrix
    get_total_nb_elements = this%n
    ! number of sub/superdiagonals is equal, every next one has 1 element less than
    ! the previous one
    do i = 1, this%kd
      get_total_nb_elements = get_total_nb_elements + 2*(this%n - i)
    end do
  end function get_total_nb_elements


  !> Checks whether the given <tt>uplo</tt> parameter is valid.
  pure logical function uplo_is_valid(uplo)
    !> uplo character to check
    character, intent(in) :: uplo

    uplo_is_valid = (uplo == "U") .or. (uplo == "L")
  end function uplo_is_valid


  !> Checks if a given position (row, col) is within the banded structure.
  !! For <tt>uplo = "U"</tt> the position is within the band if
  !! $$ max(1, col - kd) \leq row \leq col, $$
  !! for <tt>uplo = "L"</tt> the position is within the band if
  !! $$ col \leq row \leq min(n, col + kd) $$
  !! with $kd$ the number of sub/superdiagonals and $n$ the number of rows/columns.
  pure logical function is_within_band(matrix, row, col)
    !> Hermitian banded matrix structure
    class(hermitian_banded_matrix_t), intent(in) :: matrix
    !> the row index
    integer, intent(in) :: row
    !> the column index
    integer, intent(in) :: col

    if (matrix%uplo == "U") then
      is_within_band = (max(1, col - matrix%kd) <= row) .and. (row <= col)
    else  ! matrix%uplo == "L"
      is_within_band = (col <= row) .and. (row <= min(matrix%n, col + matrix%kd))
    end if
  end function is_within_band


  !> Checks if a Hermitian band matrix is compatible with another Hermitian band matrix.
  !! This implies that the following attributes should be equal:
  !!   - number of rows/columns
  !!   - number of sub/superdiagonals
  !!   - storage of upper or lower triangular part
  !!   - dimensions of the banded matrices themselves
  !! Returns <tt>.true.</tt> if all criteria are satisfied, <tt>.false.</tt> otherwise.
  pure logical function is_compatible_with(this, other)
    !> type instance
    class(hermitian_banded_matrix_t), intent(in) :: this
    !> other banded matrix
    class(hermitian_banded_matrix_t), intent(in) :: other

    is_compatible_with = ( &
      this%n == other%n &
      .and. this%kd == other%kd &
      .and. this%uplo == other%uplo &
      .and. all(shape(this%AB) == shape(other%AB)) &
    )
  end function is_compatible_with
end module mod_banded_matrix_hermitian
