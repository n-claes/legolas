module mod_solvers
  use mod_global_variables
  implicit none

  public



contains

  !> Solves the general eigenvalue problem AX = wBX with eigenvalues w
  !! using the QR-algorithm. A and B are already in tri-diagonal form.
  !! Calls zzgev from the LAPACK library.
  !! @param[in] A   Matrix A in AX = wBX
  !! @param[in] B   Matrix B in AX = wBX
  subroutine solve_QR(A, B)
    complex(dp), intent(in)  :: A(matrix_gridpts, matrix_gridpts)
    real(dp), intent(in)     :: B(matrix_gridpts, matrix_gridpts)

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

    integer                  :: i


    !! Copy contents of A and B in new arrays, as these get overwritten
    !! when calling LAPACK routines.
    A_sol = A
    B_sol = B

    !! Calculate eigenvectors or not ('N' is no, 'V' is yes)
    jobvl = 'V'
    jobvr = 'V'

    !! Get different array dimensions
    N   = matrix_gridpts
    lda = N
    ldb = N
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

    lwork = 4 * N
    allocate(work(lwork))

    call zggev(jobvl, jobvr, N, A_sol, lda, B_sol, ldb, alpha, beta, &
               vl, ldvl, vr, ldvr, work, lwork, rwork, info)

    if (info .ne. 0) then
      write(*, *) 'LAPACK routine zggev failed'
      write(*, *) 'Value for info parameter: ', info
    end if

    deallocate(vl)
    deallocate(vr)
    deallocate(work)

    ! do i = 1, matrix_gridpts
    !   write(*, *) alpha(i),  beta(i)
    ! end do

  end subroutine solve_QR

end module mod_solvers
