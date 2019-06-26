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
  !! using the QR-algorithm, the matrix B is always block-tridiagonal,
  !! symmetric and real. The problem is hence rewritten as wX = B^{-1}AX,
  !! and 'zgeev' is called from the LAPACK library.
  !! @param[in] A   Matrix A in AX = wBX
  !! @param[in] B   Matrix B in AX = wBX
  subroutine solve_QR(A, B)
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
    integer                  :: N, lda, ldb
    complex(dp)              :: A_sol(matrix_gridpts, matrix_gridpts)
    complex(dp)              :: B_sol(matrix_gridpts, matrix_gridpts)
    !! Eigenvalue variables
    complex(dp)              :: alpha(matrix_gridpts), beta(matrix_gridpts)
    !! Work variables
    integer                  :: lwork, info
    complex(dp), allocatable :: work(:)
    real(dp)                 :: rwork(8 * matrix_gridpts)

    call invert_B(B, B_inv)

    ! \TODO: code below must be reworked, call 'zgeev'



    ! !! Copy contents of A and B in new arrays, as these get overwritten
    ! !! when calling LAPACK routines.
    ! A_sol = A
    ! B_sol = B
    !
    ! !! Calculate eigenvectors or not ('N' is no, 'V' is yes)
    ! jobvl = 'V'
    ! jobvr = 'V'
    !
    ! !! Get different array dimensions
    ! N   = matrix_gridpts
    ! lda = N
    ! ldb = N
    ! if (jobvl == 'V') then
    !   ldvl = N
    ! else
    !   ldvl = 1
    ! end if
    ! if (jobvr == 'V') then
    !   ldvr = N
    ! else
    !   ldvr = 1
    ! end if
    ! allocate(vl(ldvl, N))
    ! allocate(vr(ldvr, N))
    !
    ! lwork = 4 * N
    ! allocate(work(lwork))
    !
    ! call zggev(jobvl, jobvr, N, A_sol, lda, B_sol, ldb, alpha, beta, &
    !            vl, ldvl, vr, ldvr, work, lwork, rwork, info)
    !
    ! if (info .ne. 0) then
    !   write(*, *) 'LAPACK routine zggev failed'
    !   write(*, *) 'Value for info parameter: ', info
    ! end if
    !
    ! deallocate(vl)
    ! deallocate(vr)
    ! deallocate(work)

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

end module mod_solvers
