! =============================================================================
!> Submodule containing the implementation of the QR-invert algorithm.
!! The original problem is written as a standard eigenvalue problem through
!! $$ \mathcal{B}^{-1}\mathcal{A}\textbf{X} = \omega\textbf{X}\ $$. This is
!! done using a LU decomposition via LAPACKS's <tt>zgbsv</tt>.
!! Eventually a call to LAPACK's <tt>zgeev</tt> routine is done to obtain
!! all eigenvalues and eigenvectors.
submodule (mod_solvers) smod_qr_invert
  use mod_global_variables, only: dim_subblock
  use mod_banded_matrix, only: banded_matrix_t
  use mod_transform_matrix, only: matrix_to_banded
  implicit none

contains

  !> Solves the eigenvalue problem by rewriting it to a standard form
  !! through inversion of the B-matrix.
  !! @warning Throws an error if <tt>matrix_A</tt> or <tt>matrix_B</tt>
  !!          is not a square matrix. @endwarning
  module procedure qr_invert
    !> full array containing the \(B^{-1}A\)-matrix
    complex(dp), allocatable :: array_B_invA(:, :)
    !> pivoting array
    integer, allocatable     :: ipiv(:)
    !> banded B-matrix
    type(banded_matrix_t) :: B_band
    integer     :: kd

    !> order of matrix \(B^{-1}A\)
    integer     :: N
    !> leading dimension of matrix \(B^{-1}A\)
    integer     :: ldB_invA
    !> calculate left eigenvectors if "V", omit if "N"
    character   :: jobvl
    !> leading dimension of vl
    integer     :: ldvl
    !> calculate right eigenvectors if "V", omit if "N"
    character   :: jobvr
    !> leading dimension of vr
    integer     :: ldvr
    !> dimension of work array
    integer     :: lwork
    !> work array
    complex(dp), allocatable  :: work(:)
    !> info parameter, 0 on successful exit
    integer     :: info
    !> second work array
    real(dp), allocatable     :: rwork(:)
    !> dummy for left eigenvectors, jobvl = "N" so this is never referenced
    complex(dp) :: vl(2, 2)

    kd = 2*dim_subblock+1 ! at most 2 subblocks away from diag
    call matrix_to_banded(matrix=matrix_B, subdiags=kd, superdiags=kd, banded=B_band)

    allocate(array_B_invA(matrix_A%matrix_dim, matrix_A%matrix_dim))
    call matrix_to_array(matrix=matrix_A, array=array_B_invA)
    allocate(ipiv(matrix_A%matrix_dim))
    call zgbsv( &
      B_band%n, kd, kd, N, B_band%AB, 2*kd+kd+1, ipiv, array_B_invA, ldB_invA, info &
    )
    call B_band%destroy()

    ! calculate eigenvectors, we don't use the left ones
    jobvl = "N"
    if (should_compute_eigenvectors()) then
      jobvr = "V"
    else
      jobvr = "N"
    end if
    ! set array dimensions
    N = size(array_B_invA, dim=1)
    ldB_invA = N
    ldvl = N
    ldvr = N
    ! allocate rwork array
    allocate(rwork(2 * N))
    ! get lwork
    allocate(work(1))
    call zgeev( &
      jobvl, jobvr, N, array_B_invA, ldB_invA, omega, &
      vl, ldvl, vr, ldvr, work, -1, rwork, info &
    )
    lwork = int(work(1))
    deallocate(work)
    ! allocate work array
    allocate(work(lwork))


    ! solve eigenvalue problem
    call log_message("solving evp using QR algorithm zgeev (LAPACK)", level="debug")
    call zgeev( &
      jobvl, jobvr, N, array_B_invA, ldB_invA, omega, &
      vl, ldvl, vr, ldvr, work, lwork, rwork, info &
    )
    if (info /= 0) then ! LCOV_EXCL_START
      call log_message("LAPACK routine zgeev failed!", level="warning")
      call log_message( &
        "value for the info parameter: " // str(info), &
        level="warning", &
        use_prefix=.false. &
      )
    end if ! LCOV_EXCL_STOP

    ! tear down work arrays
    deallocate(array_B_invA)
    deallocate(ipiv)
    deallocate(work)
    deallocate(rwork)
  end procedure qr_invert

end submodule smod_qr_invert
