module mod_solvers
  use mod_global_variables
  implicit none

  public



contains

  !> Solves the general eigenvalue problem AX = wBX with eigenvalues w
  !! using the QR-algorithm. A and B are already in tri-diagonal form.
  !! Calls zzgev from the LAPACK library.
  subroutine solve_QR(A, B)
    complex(dp), intent(in) :: A(matrix_gridpts, matrix_gridpts)
    real(dp), intent(in)    :: B(matrix_gridpts, matrix_gridpts)

    !! Variables for zzgev algorithm
    character               :: jobvl, jobvr
    integer                 :: n, lda, ldb, ldvl, ldvr, lwork, info, i
    complex(dp)             :: A_sol(matrix_gridpts, matrix_gridpts)
    real(dp)                :: B_sol(matrix_gridpts, matrix_gridpts)
    complex(dp)             :: alpha(matrix_gridpts), beta(matrix_gridpts)
    complex(dp)             :: vl(matrix_gridpts, matrix_gridpts)
    complex(dp)             :: vr(matrix_gridpts, matrix_gridpts)
    complex(dp)             :: work(1)
    real(dp)                :: rwork(8 * matrix_gridpts)

    !! Compute generalised eigenvectors
    jobvl = "N"
    jobvr = "N"

    N = matrix_gridpts

    !! Copy contents of A and B in new arrays, as these get overwritten
    !! when calling LAPACK routines.
    A_sol = A
    B_sol = B

    lda = N
    ldb = N
    ldvl = N
    ldvr = N
    lwork = -1

    call zggev(jobvl, jobvr, N, A_sol, lda, B_sol, ldb, alpha, beta, &
               vl, ldvl, vr, ldvr, work, lwork, rwork, info)

    if (info .ne. 0) then
      write(*, *) 'LAPACK routine zggev failed'
      write(*, *) 'Value for info parameter: ', info
    end if

    do i = 1, matrix_gridpts
      write(*, *) alpha(i),  beta(i)
    end do

  end subroutine solve_QR




end module mod_solvers
