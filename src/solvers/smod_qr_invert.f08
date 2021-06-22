! =============================================================================
!> Submodule containing the implementation of the QR-invert algorithm.
!! The original problem is written as a standard eigenvalue problem through
!! $$ \mathcal{B}^{-1}\mathcal{A}\textbf{X} = \omega\textbf{X}\ $$
!! where \(\mathcal{B}\) is always properly conditioned.
!! Various calls to <tt>mod_matrix_operations</tt> are done to invert the
!! B-matrix and to do matrix multiplications.
!! Eventually a call to LAPACK's <tt>zgeev</tt> routine is done to obtain
!! all eigenvalues and eigenvectors.
submodule (mod_solvers) smod_qr_invert
  use mod_matrix_operations, only: invert_matrix, multiply_matrices
  implicit none

contains

  !> Solves the eigenvalue problem by rewriting it to a standard form
  !! through inversion of the B-matrix.
  !! @warning Throws an error if <tt>matrix_A</tt> or <tt>matrix_B</tt>
  !!          is not a square matrix. @endwarning
  module procedure qr_invert
    !> inverse B-matrix
    real(dp)    :: B_inv(size(matrix_B, dim=1), size(matrix_B, dim=2))
    !> matrix \(B^{-1}A\)
    complex(dp) :: B_invA(size(matrix_A, dim=1), size(matrix_A, dim=2))
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

    ! do inversion of B
    call invert_matrix(matrix_B, B_inv)
    ! do matrix multiplication B^{-1}A
    call multiply_matrices(B_inv, matrix_A, B_invA)
    ! calculate eigenvectors, we don't use the left ones
    jobvl = "N"
    if (write_eigenfunctions) then
      jobvr = "V"
    else
      jobvr = "N"
    end if
    ! set array dimensions
    N = size(B_invA, dim=1)
    ldB_invA = N
    ldvl = N
    ldvr = N
    ! set work arrays
    lwork = 4 * N
    allocate(work(lwork))
    allocate(rwork(2 * N))

    ! solve eigenvalue problem
    call log_message("solving evp using QR algorithm zgeev (LAPACK)", level="debug")
    call zgeev( &
      jobvl, jobvr, N, B_invA, ldB_invA, omega, &
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

    deallocate(work)
    deallocate(rwork)

    call set_small_values_to_zero(omega)
  end procedure qr_invert

end submodule smod_qr_invert
