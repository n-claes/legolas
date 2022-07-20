! =============================================================================
!> Submodule containing the implementation of the QZ-direct algorithm.
!! We keep the general form of the eigenvalue problem
!! $$ \mathcal{A}\textbf{X} = \omega\mathcal{B}\textbf{X}\ $$
!! and solve this directly by calling LAPACK's <tt>zggev</tt> routine.
!! @note  Because the eigenvalue problem remains in general form, this
!!        routine returns the generalised eigenvectors instead of the ordinary ones.
!!        If you want eigenvectors as well, either use the QR-invert solvers
!!        or one of the ARPACK methods.
submodule (mod_solvers) smod_qz_direct
  implicit none

contains

  !> Solves the eigenvalue problem directly.
  module procedure qz_direct
    !> calculate left eigenvectors if "V", omit if "N"
    character   :: jobvl
    !> calculate right eigenvectors if "V", omit if "N"
    character   :: jobvr
    !> order of matrix_A and matrix_B
    integer     :: N
    !> leading dimension of matrix_A
    integer     :: lda
    !> full array for matrix A
    complex(dp), allocatable :: array_A(:, :)
    !> full array for matrix B
    complex(dp), allocatable :: array_B(:, :)
    !> leading dimension of (complex) matrix_B
    integer     :: ldb
    !> array for alpha, see zggem
    complex(dp) :: alpha(matrix_A%matrix_dim)
    !> array for beta, see zggem
    complex(dp) :: beta(matrix_B%matrix_dim)
    !> leading dimension of vl
    integer     :: ldvl
    !> leading dimension of vr
    integer     :: ldvr
    !> work array
    complex(dp), allocatable  :: work(:)
    !> dimension of work array
    integer     :: lwork
    !> second work array
    real(dp), allocatable     :: rwork(:)
    !> info parameter, 0 on successful exit
    integer     :: info
    !> dummy for left eigenvectors, jobvl = "N" so this is never referenced
    complex(dp) :: vl(2, 2)

    allocate(array_B(matrix_B%matrix_dim, matrix_B%matrix_dim))
    call matrix_to_array(matrix=matrix_B, array=array_B)
    allocate(array_A(matrix_A%matrix_dim, matrix_A%matrix_dim))
    call matrix_to_array(matrix=matrix_A, array=array_A)

    !> @warning The LAPACK routine <tt>zggev</tt> returns the _generalised_
    !! eigenvectors, which are different from the ordinary ones returned
    !! by the QR algorithm <tt>zgeev</tt>.
    !! For the moment we don't calculate eigenvectors with the direct solver until
    !! we get consistent results between QZ and QR. In the meantime,
    !! use <tt>solver = "QR-invert"</tt> if eigenvectors are needed. @endwarning
    jobvl = "N"
    jobvr = "N"
    if (write_eigenfunctions) then  ! LCOV_EXCL_START
      call log_message( &
        "eigenvector calculations with the direct QZ solver are disabled for now,", &
        level="warning" &
      )
      call log_message( &
        "use the QR-invert solver instead", &
        level="warning", &
        use_prefix=.false. &
      )
    end if  ! LCOV_EXCL_STOP
    ! set array dimensions
    N = matrix_A%matrix_dim
    lda = N
    ldb = N
    ldvl = N
    ldvr = N
    ! set work arrays
    lwork = 4 * N
    allocate(work(lwork))
    allocate(rwork(8 * N))

    ! solve eigenvalue problem
    call log_message("solving evp using QZ algorithm zggev (LAPACK)", level="debug")
    call zggev( &
      jobvl, jobvr, N, array_A, lda, array_B, ldb, &
      alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info &
    )
    if (info /= 0) then ! LCOV_EXCL_START
      call log_message("LAPACK routine zggev failed!", level="warning")
      call log_message( &
        "value for the info parameter: " // str(info), &
        level="warning", &
        use_prefix=.false. &
      )
    end if ! LCOV_EXCL_STOP
    omega = alpha / beta

    deallocate(work)
    deallocate(rwork)
    deallocate(array_B)
    deallocate(array_A)

    call set_small_values_to_zero(omega)
  end procedure qz_direct
end submodule smod_qz_direct
