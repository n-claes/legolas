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

  private

  public  :: solve_QR

contains

  !> Solves the general eigenvalue problem AX = wBX with eigenvalues w
  !! using the QR-algorithm. The problem is rewritten as wX = B^{-1}AX,
  !! and 'zgeev' is called from the LAPACK library.
  !! @param[in] A   Matrix A in AX = wBX
  !! @param[in] B   Matrix B in AX = wBX
  !! @param[out] omega  Array of size matrix_gridpts containing the eigenvalues
  subroutine solve_QR(A, B, omega)
    complex(dp), intent(in)  :: A(matrix_gridpts, matrix_gridpts)
    real(dp), intent(in)     :: B(matrix_gridpts, matrix_gridpts)
    real(dp)                 :: B_inv(matrix_gridpts, matrix_gridpts)

    !! Left eigenvector variables
    character                :: jobvl
    complex(dp), allocatable :: vl(:, :)
    integer                  :: ldvl
    !! Right eigenvector variables
    character                :: jobvr
    complex(dp), allocatable :: vr(:, :)
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
    jobvl = 'V'
    jobvr = 'V'

    !! Array dimensions
    N       = matrix_gridpts
    ldB_invA = N
    if (jobvl == 'V') then
      ldvl = N
    else
      ldvl = 1
    end if
    if (jobvr == 'V') then
      ldvr = N
    else
      ldvr = 1
    end if

    allocate(vl(ldvl, N))
    allocate(vr(ldvr, N))

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

    deallocate(vl)
    deallocate(vr)
    deallocate(work)
    deallocate(rwork)

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


    K   = matrix_gridpts
    ldB_inv  = K
    ldA      = K
    ldB_invA = K

    alpha = (1.0d0, 0.0d0)
    beta  = (0.0d0, 0.0d0)

    !! 'zgemm' performs one of the matrix-matrix operations
    !! C := alpha*op(A)*op(B) + beta*C
    !! In this case, alpha = 1, beta = 0, op = 'N' (so no transp. or conj.)
    call zgemm('N', 'N', K, K, K, alpha, B_inv, ldB_inv, A, ldA, &
               beta, B_invA, ldB_invA)

  end subroutine get_B_invA

end module mod_solvers
