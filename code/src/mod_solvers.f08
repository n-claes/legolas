!
! MODULE: mod_solvers
!
!> @author
!> Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Module containing methods to solve the eigenvalue problem AX = wBX.
!! Does calls to the BLAS and LAPACK libraries.
!
module mod_solvers
  use mod_global_variables
  implicit none

  public

contains

  !> Solves the general eigenvalue problem AX = wBX with eigenvalues w
  !! using the QR-algorithm. The problem is rewritten as wX = B^{-1}AX,
  !! and 'zgeev' is called from the LAPACK library.
  !! @param[in] A   Matrix A in AX = wBX
  !! @param[in] B   Matrix B in AX = wBX
  !! @param[out] omega  Array of size matrix_gridpts containing the eigenvalues
  !! @param[out] vl The left eigenvectors, only calculated when
  !!                write_eigenfunctions is .true., otherwise zero.
  !! @param[out] vr The right eigenvectors, only calculated when
  !!                write_eigenfunctions is .true., otherwise zero.
  subroutine solve_QR(A, B, omega, vl, vr)
    use mod_check_values

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
    !call get_B_invA_matmul(B_inv, A, B_invA)
    call get_B_invA(B_inv, A, B_invA)

    !! Calculate eigenvectors or not ('N' is no, 'V' is yes)
    if (write_eigenvectors) then
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

    if (info .ne. 0) then
      write(*, *) 'LAPACK routine zggev failed'
      write(*, *) 'Value for info parameter: ', info
    end if

    if (.not. write_eigenvectors) then
      vl = (0.0d0, 0.0d0)
      vr = (0.0d0, 0.0d0)
    end if

    deallocate(work)
    deallocate(rwork)

    call check_small_values(omega)
    if (write_eigenvectors) then
      call check_small_values_matrix(vl)
      call check_small_values_matrix(vr)
    end if

  end subroutine solve_QR


  !> Inverts the matrix B using LAPACK routines. First a LU-factorisation is
  !! done by 'dgetrf', which is then used by 'dgetri' to calculate the inverse.
  !! @param[in] B       Matrix B from AX = wBX.
  !!                    Should be real, symmetric and block-tridiagonal.
  !! @param[out] B_inv  Inverse of B
  subroutine invert_B(B, B_inv)
    real(dp), intent(in)  :: B(matrix_gridpts, matrix_gridpts)
    real(dp), intent(out) :: B_inv(matrix_gridpts, matrix_gridpts)
    integer               :: N, ldb, lwork, info
    integer, allocatable  :: ipiv(:)
    real(dp), allocatable :: work(:)

    !! Copy B into B_inv
    B_inv = B

    N   = matrix_gridpts
    ldb = N

    lwork = 4*N

    allocate(ipiv(N))
    allocate(work(lwork))

    ! Calculate pivot indices
    call dgetrf(N, N, B_inv, ldb, ipiv, info)
    if (info .ne. 0) then
      write(*, *) 'LU factorisation of matrix B failed'
      write(*, *) 'Value for info parameter: ', info
    end if

    call dgetri(N, B_inv, ldb, ipiv, work, lwork, info)

    if (info .ne. 0) then
      write(*, *) 'Inversion of matrix B failed'
      write(*, *) 'Value for info parameter: ', info
    end if

    deallocate(ipiv)
    deallocate(work)

  end subroutine invert_B


  !> Does the matrix multiplication B^{-1}A. The LAPACK routine 'zgemm' is
  !! used instead of the build-in 'matmul' function for speedup and efficiency.
  !! @param[in]   B_inv   The inverse of the matrix B
  !! @param[in]   A       The matrix A
  !! @param[out]  B_invA  The result of the matrix multiplication B_inv * A
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

    call zgemm('N', 'N', K, K, K, alpha, B_inv_cplx, ldB_inv, A, ldA, &
               beta, B_invA, ldB_invA)

  end subroutine get_B_invA


  !> Also does the matrix multiplication B^{-1}A, but using Fortran's
  !! build-in 'matmul' function. Used for comparison with the case above.
  !! @param[in]   B_inv   The inverse of the matrix B
  !! @param[in]   A       The matrix A
  !! @param[out]  B_invA  The result of the matrix multiplication B_inv * A
  subroutine get_B_invA_matmul(B_inv, A, B_invA)
    real(dp), intent(in)      :: B_inv(matrix_gridpts, matrix_gridpts)
    complex(dp), intent(in)   :: A(matrix_gridpts, matrix_gridpts)
    complex(dp), intent(out)  :: B_invA(matrix_gridpts, matrix_gridpts)

    B_invA = matmul(B_inv, A)

  end subroutine get_B_invA_matmul


  !> Also does the matrix multiplication B^{-1}A, but is manually implemented
  !! using three do-loops. Used for comparison with the other methods.
  !! @param[in]   B_inv   The inverse of the matrix B
  !! @param[in]   A       The matrix A
  !! @param[out]  B_invA  The result of the matrix multiplication B_inv * A
  subroutine get_B_invA_manual(B_inv, A, B_invA)
    real(dp), intent(in)      :: B_inv(matrix_gridpts, matrix_gridpts)
    complex(dp), intent(in)   :: A(matrix_gridpts, matrix_gridpts)
    complex(dp), intent(out)  :: B_invA(matrix_gridpts, matrix_gridpts)

    integer                   :: i, j, k

    B_invA = (0.0d0, 0.0d0)

    do i = 1, matrix_gridpts
      do j = 1, matrix_gridpts
        do k = 1, matrix_gridpts
          B_invA(i, j) = B_invA(i, j) + B_inv(i, k) * A(k, j)
        end do
      end do
    end do

  end subroutine get_B_invA_manual

end module mod_solvers
