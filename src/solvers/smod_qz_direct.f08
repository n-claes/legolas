! =============================================================================
!> Submodule containing the implementation of the QZ-direct algorithm.
!! We keep the general form of the eigenvalue problem
!! $$ \mathcal{A}\textbf{X} = \omega\mathcal{B}\textbf{X}\ $$
!! and solve this directly by calling LAPACK's <tt>zggev3</tt> routine.
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

    jobvl = "N"
    jobvr = "N"
    if (settings%io%should_compute_eigenvectors()) jobvr = "V"
    ! set array dimensions
    N = matrix_A%matrix_dim
    lda = N
    ldb = N
    ldvl = N
    ldvr = N
    ! allocate rwork array
    allocate(rwork(8 * N))
    ! get lwork
    allocate(work(1))
    call zggev( &
      jobvl, jobvr, N, array_A, lda, array_B, ldb, &
      alpha, beta, vl, ldvl, vr, ldvr, work, -1, rwork, info &
    )
    lwork = int(work(1))
    deallocate(work)
    ! allocate work array
    allocate(work(lwork))

    ! solve eigenvalue problem
    call logger%debug("solving evp using QZ algorithm zggev (LAPACK)")
    call zggev( &
      jobvl, jobvr, N, array_A, lda, array_B, ldb, &
      alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info &
    )
    if (info /= 0) then ! LCOV_EXCL_START
      call logger%warning("LAPACK routine zggev failed!")
      call logger%warning("value for the info parameter: " // str(info))
    end if ! LCOV_EXCL_STOP
    omega = alpha / beta

    deallocate(work)
    deallocate(rwork)
    deallocate(array_B)
    deallocate(array_A)

    call set_small_values_to_zero(omega)
  end procedure qz_direct
end submodule smod_qz_direct
