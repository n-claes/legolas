! =============================================================================
!> @brief   Module containing the eigenvalue problem solvers.
!! @details This module handles everything related to solving the eigenvalue problem
!!          \f$\mathcal{A}\textbf{X} = \omega\mathcal{B}\textbf{X}\f$,
!!          that is, matrix inversion, matrix multiplication and the QR-algorithm.
!!          Calls to LAPACK and BLAS are done in this module.
module mod_solvers
  use mod_global_variables, only: dp, matrix_gridpts
  use mod_logging, only: log_message, char_log, int_fmt
  implicit none

  private

  public :: solve_QR
  public :: invert_B
  public :: get_B_invA

contains


  !> @brief   QR-solver for the eigenvalue problem.
  !! @details Solves a complex non-Hermitian general eigenvalue problem using
  !!          the QR-algorithm. The matrix B is inverted in doing so.
  !!          Also calculates the eigenvectors if requested.
  !!          Uses the \p zgeev routine from LAPACK.
  !! @exception Warning If the \p zgeev routine fails.
  !! @note  All eigenvalues below \p DP_LIMIT are set to zero. The real and imaginary
  !!        parts are handled separately.
  !! @note  The eigenvectors are only calculated if \p write_eigenfunctions is \p True.
  !!        Otherwise an array of zeroes is returned.
  !! @param[in] A     matrix A
  !! @param[in] B     matrix B
  !! @param[out] omega  array containing the calculated eigenvalues
  !! @param[out] vl     array containing the left eigenvectors
  !! @param[out] vr     array containing the right eigenvectors
  subroutine solve_QR(A, B, omega, vl, vr)
    use mod_global_variables, only: write_eigenfunctions
    use mod_check_values, only: check_small_values
    use mod_logging, only: log_message, int_fmt, char_log

    complex(dp), intent(in)  :: A(matrix_gridpts, matrix_gridpts)
    real(dp), intent(in)     :: B(matrix_gridpts, matrix_gridpts)
    real(dp)                 :: B_inv(matrix_gridpts, matrix_gridpts)

    !! Left eigenvector variables
    character                :: jobvl
    complex(dp), intent(out) :: vl(matrix_gridpts, matrix_gridpts)
    integer                  :: ldvl
    !! Right eigenvector variables
    character                :: jobvr
    complex(dp), intent(out) :: vr(matrix_gridpts, matrix_gridpts)
    integer                  :: ldvr
    !! Matrix variables
    integer                  :: N, ldB_invA
    complex(dp)              :: B_invA(matrix_gridpts, matrix_gridpts)
    !! Eigenvalue variables
    complex(dp), intent(out) :: omega(matrix_gridpts)
    !! Work variables
    integer                  :: lwork, info
    complex(dp), allocatable :: work(:)
    real(dp), allocatable    :: rwork(:)

    !! Invert matrix B
    call invert_B(B, B_inv)
    !! Matrix multiplication B^{-1} * A
    call get_B_invA(B_inv, A, B_invA)

    !! Calculate eigenvectors or not ('N' is no, 'V' is yes)
    if (write_eigenfunctions) then
      jobvl = 'V'
      jobvr = 'V'
    else
      jobvl = 'N'
      jobvr = 'N'
    end if

    !! Array dimensions
    N       = matrix_gridpts
    ldB_invA = N
    ldvl    = N
    ldvr    = N

    !! Size or work array
    lwork = 4*N
    allocate(work(lwork))
    allocate(rwork(2*N))

    call zgeev(jobvl, jobvr, N, B_invA, ldB_invA, omega, vl, ldvl, &
               vr, ldvr, work, lwork, rwork, info)

    if (info /= 0) then
      write(char_log, int_fmt) info
      call log_message('LAPACK routine zgeev failed!', level='warning')
      call log_message('value for the info parameter: ' // adjustl(char_log), level='warning')
    end if

    deallocate(work)
    deallocate(rwork)

    call check_small_values(omega)
  end subroutine solve_QR


  !> @brief   Handles the inversion of the B-matrix.
  !! @details Inverts the B-matrix using LAPACK routines. First a LU-factorisation is performed using
  !!          \p dgetrf, after which inversion is done through the routine \p dgetri.
  !! @exception Warning   If the \p dgetri routine fails.
  !! @exception Warning   If the \p dgetrf routine fails.
  !! @param[in]  B      the B-matrix
  !! @param[out] B_inv  the inverted B-matrix
  subroutine invert_B(B, B_inv)
    real(dp), intent(in)  :: B(matrix_gridpts, matrix_gridpts)
    real(dp), intent(out) :: B_inv(matrix_gridpts, matrix_gridpts)
    integer               :: N, ldb, lwork, info
    integer, allocatable  :: ipiv(:)
    real(dp), allocatable :: work(:)

    call log_message("inverting B-matrix", level='debug')

    !! Copy B into B_inv
    B_inv = B

    N   = matrix_gridpts
    ldb = N
    lwork = 4*N

    allocate(ipiv(N))
    allocate(work(lwork))

    ! Calculate pivot indices
    call log_message("LU factorisation of B using dgetrf", level='debug')
    call dgetrf(N, N, B_inv, ldb, ipiv, info)
    if (info /= 0) then
      write(char_log, int_fmt) info
      call log_message("LU factorisation of B failed. Value info: " // trim(char_log), level='warning')
    end if

    call log_message("inverting B using dgetri", level='debug')
    call dgetri(N, B_inv, ldb, ipiv, work, lwork, info)

    if (info /= 0) then
      write(char_log, int_fmt) info
      call log_message("inversion of B failed. Value info: " // adjustl(char_log), level='warning')
    end if

    deallocate(ipiv)
    deallocate(work)
  end subroutine invert_B


  !> @brief   Handles matrix multiplication of B inverse with A.
  !! @details Does the matrix multiplication \f$\mathcal{B}^{-1}\mathcal{A}\f$ using LAPACK.
  !!          The routine \p zgemm is used instead of the Fortran builtin \p matmul for efficiency.
  !! @param[in]   B_inv   the inverse of B
  !! @param[in]   A       the matrix A
  !! @param[out]  B_invA  the result of the matrix multiplication <tt>B_inv*A</tt>
  subroutine get_B_invA(B_inv, A, B_invA)
    real(dp), intent(in)      :: B_inv(matrix_gridpts, matrix_gridpts)
    complex(dp), intent(in)   :: A(matrix_gridpts, matrix_gridpts)
    complex(dp), intent(out)  :: B_invA(matrix_gridpts, matrix_gridpts)

    integer                   :: K, ldB_inv, ldA, ldB_invA
    complex(dp)               :: alpha, beta
    complex(dp)               :: B_inv_cplx(matrix_gridpts, matrix_gridpts)

    !! 'zgemm' performs one of the matrix-matrix operations
    !! C := alpha*op(A)*op(B) + beta*C
    !! In this case, alpha = 1, beta = 0, op = 'N' (so no transp. or conj.)
    K   = matrix_gridpts
    ldB_inv  = K
    ldA      = K
    ldB_invA = K

    alpha = (1.0d0, 0.0d0)
    beta  = (0.0d0, 0.0d0)

    !! The input matrix B HAS TO BE COMPLEX (in our case it is real).
    !! Hence convert it to a complex matrix, otherwise results are wrong.
    B_inv_cplx = B_inv * (1.0d0, 0.0d0)

    call log_message("multiplying B_inv * A", level='debug')
    call zgemm('N', 'N', K, K, K, alpha, B_inv_cplx, ldB_inv, A, ldA, &
               beta, B_invA, ldB_invA)
  end subroutine get_B_invA

end module mod_solvers
