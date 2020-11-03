module mod_matrix_operations
  use mod_global_variables, only: dp
  use mod_check_values, only: matrix_is_square
  use mod_logging, only: log_message, char_log, int_fmt
  implicit none

  !> integer used to set matrix dimensions
  integer   :: N
  !> integer used to set leading dimension of matrix
  integer   :: ldm1
  !> integer used to set leading dimension of (second) matrix
  integer   :: ldm2
  !> integer used to set size of work array
  integer   :: lwork
  !> integer used for successful exits
  integer   :: info
  !> array for pivot indices
  integer, allocatable  :: ipiv(:)
  !> array for work
  real(dp), allocatable :: work(:)

  private

  !> Interface to invert a matrix
  interface invert_matrix
    module procedure invert_matrix_real
  end interface invert_matrix

  !> Interface to do matrix multiplication
  interface multiply_matrices
    module procedure rmat_x_cmat
    module procedure cmat_x_carr
  end interface multiply_matrices

  public  :: invert_matrix
  public  :: multiply_matrices

contains


  !> Handles the inversion of a real square matrix using LAPACK routines.
  !! First a LU-factorisation is performed using <tt>dgetrf</tt>, then inversion is
  !! done using <tt>dgetri</tt>.
  !! @warning Throws a warning if <tt>dgetrf</tt> or <tt>dgetri</tt> fails. @endwarning
  !! @warning Throws an error if the matrix is not square.
  subroutine invert_matrix_real(mat, mat_inv)
    !> real matrix to invert
    real(dp), intent(in)  :: mat(:, :)
    !> inverted matrix
    real(dp), intent(out) :: mat_inv(:, :)

    if (.not. matrix_is_square(mat)) then
      call log_message("trying to invert but matrix is not square!", level="error")
    end if

    mat_inv = mat

    N = size(mat, dim=1)
    ldm1 = N
    lwork = 4 * N
    allocate(ipiv(N))
    allocate(work(lwork))

    ! calculate pivot indices
    call log_message("LU factorisation of matrix using dgetrf", level="debug")
    call dgetrf(N, N, mat_inv, ldm1, ipiv, info)
    if (info /= 0) then
      write(char_log, int_fmt) info
      call log_message( &
        "LU factorisation of matrix failed. Value info: " // adjustl(char_log), &
        level="warning" &
      )
    end if
    ! invert matrix
    call log_message("inverting matrix using dgetri", level="debug")
    call dgetri(N, mat_inv, ldm1, ipiv, work, lwork, info)
    if (info /= 0) then
      write(char_log, int_fmt) info
      call log_message( &
        "inversion of matrix failed. Value info: " // adjustl(char_log), &
        level="warning" &
      )
    end if

    deallocate(ipiv)
    deallocate(work)
  end subroutine invert_matrix_real


  !> Matrix multiplication using LAPACK routines, multiplies
  !! a real with a complex matrix.
  subroutine rmat_x_cmat(mat1, mat2, mat_out)
    !> first matrix (left side)
    real(dp), intent(in)      :: mat1(:, :)
    !> second matrix (right side)
    complex(dp), intent(in)   :: mat2(:, :)
    !> matrix multiplication mat1 x mat2
    complex(dp), intent(out)  :: mat_out(:, :)

    !> number of rows in mat1
    integer   :: rows_mat1
    !> number of columns in mat1
    integer   :: cols_mat1
    !> number of rows in mat2
    integer   :: rows_mat2
    !> number of columns in mat2
    integer   :: cols_mat2
    !> complex analog of mat1
    complex(dp) :: mat1_c(size(mat1, dim=1), size(mat1, dim=2))
    !> scalar alpha, see zgemm
    complex(dp) :: alpha
    !> scalar beta, see zgemm
    complex(dp) :: beta

    ldm1 = size(mat1, dim=1)
    ldm2 = size(mat2, dim=1)
    rows_mat1 = size(mat1, dim=1)
    cols_mat1 = size(mat1, dim=2)
    rows_mat2 = size(mat2, dim=1)
    cols_mat2 = size(mat2, dim=2)

    ! check compatibility
    if (cols_mat1 /= rows_mat2) then
      call log_message( &
        "incompatible matrices during matrix multiplication", &
        level="error" &
      )
    end if

    !> @note <tt>zgemm</tt> performs one of the matrix-matrix operations
    !! $$ C := \alpha*op(A)*op(B) + \beta*C $$
    !! In this case, \(\alpha = 1\), \(\beta = 0\), op = 'N'
    !! (so no transpose or conjugate).
    alpha = (1.0d0, 0.0d0)
    beta  = (0.0d0, 0.0d0)
    ! make mat1 complex
    mat1_c = mat1 * (1.0d0, 0.0d0)
    ! do multiplication
    call zgemm( &
      "N", "N", rows_mat1, cols_mat2, rows_mat2, &
      alpha, mat1_c, ldm1, mat2, ldm2, beta, mat_out, ldm1 &
    )
  end subroutine rmat_x_cmat


  !> Matrix multiplication using LAPACK routines,
  !! multiplies a complex matrix with a complex vector
  subroutine cmat_x_carr(mat, vec, vec_out)
    !> matrix (left side)
    complex(dp), intent(in)   :: mat(:, :)
    !> column vector (right side)
    complex(dp), intent(in)   :: vec(:)
    !> matrix multiplication mat x vec, yields vector
    complex(dp), intent(out)  :: vec_out(:)

    !> number of rows in mat
    integer   :: rows_mat
    !> number of cols in mat
    integer   :: cols_mat
    !> number of rows in vec
    integer   :: rows_vec
    !> number of cols in vec
    integer   :: cols_vec
    !> scalar alpha, see zgemm
    complex(dp) :: alpha
    !> scalar beta, see zgemm
    complex(dp) :: beta

    ldm1 = size(mat, dim=1)
    ldm2 = size(vec, dim=1)
    rows_mat = size(mat, dim=1)
    cols_mat = size(mat, dim=2)
    rows_vec = size(vec, dim=1)
    cols_vec = 1

    ! check compatibility
    if (cols_mat /= rows_vec) then
      call log_message( &
        "incompatible matrices during matrix multiplication", &
        level="error" &
      )
    end if

    !> @note <tt>zgemm</tt> performs one of the matrix-matrix operations
    !! $$ C := \alpha*op(A)*op(B) + \beta*C $$
    !! In this case, \(\alpha = 1\), \(\beta = 0\), op = 'N'
    !! (so no transpose or conjugate).
    alpha = (1.0d0, 0.0d0)
    beta  = (0.0d0, 0.0d0)
    ! do multiplication
    call zgemm( &
      "N", "N", rows_mat, cols_vec, cols_mat, &
      alpha, mat, ldm1, vec, ldm2, beta, vec_out, ldm1 &
    )
  end subroutine cmat_x_carr

end module mod_matrix_operations