! =============================================================================
!> Module containing various matrix operations (inverting, multiplying, etc).
!! Does calls to relevant BLAS or LAPACK routines.
module mod_matrix_operations
  use mod_global_variables, only: dp
  use mod_logging, only: log_message, str
  implicit none

  !> integer used to set leading dimension of matrix
  integer   :: rows_mat1
  !> integer used to set trailing dimension of matrix
  integer   :: cols_mat1
  !> integer used to set leading dimension of second matrix
  integer   :: rows_mat2
  !> integer used to set trailing simension of second matrix
  integer   :: cols_mat2
  !> integer used to set size of work array
  integer   :: lwork
  !> integer used for successful exits
  integer   :: info

  private

  !> Interface to invert a matrix
  interface invert_matrix
    module procedure invert_matrix_real
    module procedure invert_matrix_complex
  end interface invert_matrix

  !> Interface to do matrix multiplication
  interface multiply_matrices
    module procedure rmat_x_cmat
    module procedure rmat_x_cvec
    module procedure cmat_x_rmat
    module procedure cmat_x_cvec
  end interface multiply_matrices

  public  :: invert_matrix
  public  :: multiply_matrices

contains


  !> Handles the inversion of a real square matrix using LAPACK routines.
  !! First a LU-factorisation is performed using <tt>dgetrf</tt>, then inversion is
  !! done using <tt>dgetri</tt>.
  !! @warning Throws a warning if <tt>dgetrf</tt> or <tt>dgetri</tt> fails. @endwarning
  !! @warning Throws an error if the matrix is not square. @endwarning
  subroutine invert_matrix_real(mat, mat_inv)
    !> real matrix to invert
    real(dp), intent(in)  :: mat(:, :)
    !> inverted matrix
    real(dp), intent(out) :: mat_inv(:, :)
    !> array for pivot indices
    integer, allocatable  :: ipiv(:)
    !> array for work
    real(dp), allocatable :: work(:)

    if (size(mat, dim=1) /= size(mat, dim=2)) then
      call log_message("trying to invert but matrix is not square!", level="error")
      return
    end if

    mat_inv = mat
    rows_mat1 = size(mat, dim=1)
    lwork = 4 * rows_mat1
    allocate(ipiv(rows_mat1))
    allocate(work(lwork))
    ! calculate pivot indices
    call log_message("LU factorisation of matrix using dgetrf", level="debug")
    call dgetrf(rows_mat1, rows_mat1, mat_inv, rows_mat1, ipiv, info)
    if (info /= 0) then
      call log_message( &
        "LU factorisation of matrix failed. Value info: " // str(info), &
        level="warning" &
      )
    end if
    ! invert matrix
    call log_message("inverting matrix using dgetri", level="debug")
    call dgetri(rows_mat1, mat_inv, rows_mat1, ipiv, work, lwork, info)
    if (info /= 0) then
      call log_message( &
        "inversion of matrix failed. Value info: " // str(info), &
        level="warning" &
      )
    end if
    deallocate(ipiv)
    deallocate(work)
  end subroutine invert_matrix_real


  !> Handles the inversion of a complex square matrix using LAPACK routines.
  !! First a LU-factorisation is performed using <tt>zgetrf</tt>, then inversion is
  !! done using <tt>zgetri</tt>.
  !! @warning Throws a warning if <tt>zgetrf</tt> or <tt>zgetri</tt> fails. @endwarning
  !! @warning Throws an error if the matrix is not square. @endwarning
  subroutine invert_matrix_complex(mat, mat_inv)
    !> complex matrix to invert
    complex(dp), intent(in)  :: mat(:, :)
    !> inverted matrix
    complex(dp), intent(out) :: mat_inv(:, :)
    !> array for pivot indices
    integer, allocatable     :: ipiv(:)
    !> array for work
    complex(dp), allocatable :: work(:)

    if (size(mat, dim=1) /= size(mat, dim=2)) then
      call log_message("trying to invert but matrix is not square!", level="error")
      return
    end if
    mat_inv = mat
    rows_mat1 = size(mat, dim=1)
    lwork = 4 * rows_mat1
    allocate(ipiv(rows_mat1))
    allocate(work(lwork))
    ! calculate pivot indices
    call log_message("LU factorisation of matrix using zgetrf", level="debug")
    call zgetrf(rows_mat1, rows_mat1, mat_inv, rows_mat1, ipiv, info)
    if (info /= 0) then
      call log_message( &
        "LU factorisation of matrix failed. Value info: " // str(info), &
        level="warning" &
        )
    end if
    ! invert matrix
    call log_message("inverting matrix using zgetri", level="debug")
    call zgetri(rows_mat1, mat_inv, rows_mat1, ipiv, work, lwork, info)
    if (info /= 0) then
      call log_message( &
        "inversion of matrix failed. Value info: " // str(info), &
        level="warning" &
        )
    end if
    deallocate(ipiv)
    deallocate(work)
  end subroutine invert_matrix_complex


  !> Matrix multiplication using the LAPACK routine <tt>zgemm</tt>, multiplies
  !! a real with a complex matrix.
  subroutine rmat_x_cmat(mat1, mat2, mat_out)
    !> first matrix (left side)
    real(dp), intent(in)      :: mat1(:, :)
    !> second matrix (right side)
    complex(dp), intent(in)   :: mat2(:, :)
    !> matrix multiplication mat1 x mat2
    complex(dp), intent(out)  :: mat_out(:, :)
    !> scalar alpha, see zgemm
    complex(dp) :: alpha
    !> scalar beta, see zgemm
    complex(dp) :: beta

    rows_mat1 = size(mat1, dim=1)
    cols_mat1 = size(mat1, dim=2)
    rows_mat2 = size(mat2, dim=1)
    cols_mat2 = size(mat2, dim=2)
    call check_matrix_compatibility(cols_mat1, rows_mat2)
    alpha = (1.0d0, 0.0d0)
    beta  = (0.0d0, 0.0d0)
    ! do multiplication, zgemm needs two complex matrices
    call zgemm( &
      "N", "N", rows_mat1, cols_mat2, rows_mat2, alpha, mat1 * (1.0d0, 0.0d0), &
      rows_mat1, mat2, rows_mat2, beta, mat_out, rows_mat1 &
    )
  end subroutine rmat_x_cmat


  !> Matrix multiplication using the LAPACK routine <tt>zgemm</tt>,
  !! multiplies a real matrix with a complex vector
  subroutine rmat_x_cvec(mat, vec, vec_out)
    !> matrix (left side)
    real(dp), intent(in)      :: mat(:, :)
    !> column vector (right side)
    complex(dp), intent(in)   :: vec(:)
    !> matrix multiplication mat x vec, yields vector
    complex(dp), intent(out)  :: vec_out(:)
    !> scalar alpha, see zgemm
    complex(dp) :: alpha
    !> scalar beta, see zgemm
    complex(dp) :: beta

    rows_mat1 = size(mat, dim=1)
    cols_mat1 = size(mat, dim=2)
    rows_mat2 = size(vec, dim=1)
    cols_mat2 = 1
    call check_matrix_compatibility(cols_mat1, rows_mat2)
    alpha = (1.0d0, 0.0d0)
    beta  = (0.0d0, 0.0d0)
    call zgemm( &
      "N", "N", rows_mat1, cols_mat2, cols_mat1, alpha, mat * (1.0d0, 0.0d0), &
      rows_mat1, vec, rows_mat2, beta, vec_out, rows_mat1 &
    )
  end subroutine rmat_x_cvec


  !> Matrix multiplication using the LAPACK routine <tt>zgemm</tt>,
  !! multiplies a complex with a real matrix.
  subroutine cmat_x_rmat(mat1, mat2, mat_out)
    !> first matrix (left side)
    complex(dp), intent(in)   :: mat1(:, :)
    !> second matrix (right side)
    real(dp), intent(in)      :: mat2(:, :)
    !> matrix multiplication mat1 x mat2
    complex(dp), intent(out)  :: mat_out(:, :)
    !> scalar alpha, see zgemm
    complex(dp) :: alpha
    !> scalar beta, see zgemm
    complex(dp) :: beta

    rows_mat1 = size(mat1, dim=1)
    cols_mat1 = size(mat1, dim=2)
    rows_mat2 = size(mat2, dim=1)
    cols_mat2 = size(mat2, dim=2)
    call check_matrix_compatibility(cols_mat1, rows_mat2)
    alpha = (1.0d0, 0.0d0)
    beta  = (0.0d0, 0.0d0)
    ! do multiplication, zgemm needs two complex matrices
    call zgemm( &
      "N", "N", rows_mat1, cols_mat2, rows_mat2, alpha, mat1, &
      rows_mat1, mat2 * (1.0d0, 0.0d0), rows_mat2, beta, mat_out, rows_mat1 &
    )
  end subroutine cmat_x_rmat


  !> Matrix multiplication using the LAPACK routine <tt>zgemm</tt>,
  !! multiplies a complex matrix with a complex vector
  subroutine cmat_x_cvec(mat, vec, vec_out)
    !> matrix (left side)
    complex(dp), intent(in)   :: mat(:, :)
    !> column vector (right side)
    complex(dp), intent(in)   :: vec(:)
    !> matrix multiplication mat x vec, yields vector
    complex(dp), intent(out)  :: vec_out(:)
    !> scalar alpha, see zgemm
    complex(dp) :: alpha
    !> scalar beta, see zgemm
    complex(dp) :: beta

    rows_mat1 = size(mat, dim=1)
    cols_mat1 = size(mat, dim=2)
    rows_mat2 = size(vec, dim=1)
    cols_mat2 = 1
    call check_matrix_compatibility(cols_mat1, rows_mat2)
    alpha = (1.0d0, 0.0d0)
    beta  = (0.0d0, 0.0d0)
    ! do multiplication
    call zgemm( &
      "N", "N", rows_mat1, cols_mat2, cols_mat1, &
      alpha, mat, rows_mat1, vec, rows_mat2, beta, vec_out, rows_mat1 &
    )
  end subroutine cmat_x_cvec


  !> Checks if two matrices (or a matrix and a vector) are compatible
  !! for matrix multiplication.
  subroutine check_matrix_compatibility(cols_mat1, rows_mat2)
    !> the number of columns of the first matrix
    integer, intent(in) :: cols_mat1
    !> the number of rows of the second matrix
    integer, intent(in) :: rows_mat2

    if (cols_mat1 /= rows_mat2) then
      call log_message( &
        "incompatible matrix multiplication: (. x " // str(cols_mat1) // ") x (" &
        // str(rows_mat2) // " x .)", &
        level="error" &
      )
    end if
  end subroutine check_matrix_compatibility

end module mod_matrix_operations