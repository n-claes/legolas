submodule (mod_solvers) smod_qr_invert
  use mod_matrix_operations, only: invert_matrix, multiply_matrices
  implicit none

contains

  module subroutine qr_invert(matrix_A, matrix_B, omega, vl, vr)
    !> matrix A
    complex(dp), intent(in)   :: matrix_A(:, :)
    !> matrix B
    real(dp), intent(in)      :: matrix_B(:, :)
    !> array with calculated eigenvalues
    complex(dp), intent(out)  :: omega(:)
    !> array with left eigenvectors
    complex(dp), intent(out)  :: vl(:, :)
    !> array with right eigenvectors
    complex(dp), intent(out)  :: vr(:, :)

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

    call log_message("using QR-invert", level="debug")

    ! check that we have square matrices
    if (size(matrix_A, dim=1) /= size(matrix_A, dim=2)) then
      call log_message("qr_invert: A matrix is not square!", level="error")
      return
    end if
    if (size(matrix_B, dim=1) /= size(matrix_B, dim=2)) then
      call log_message("qr_invert: B matrix is not square!", level="error")
      return
    end if

    ! do inversion of B
    call invert_matrix(matrix_B, B_inv)
    ! do matrix multiplication B^{-1}A
    call multiply_matrices(B_inv, matrix_A, B_invA)
    ! calculate eigenvectors, we don't use the right ones
    jobvr = "N"
    if (write_eigenfunctions) then
      jobvl = "V"
    else
      jobvl = "N"
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
    call zgeev( &
      jobvl, jobvr, N, B_invA, ldB_invA, omega, &
      vl, ldvl, vr, ldvr, work, lwork, rwork, info &
    )
    if (info /= 0) then
      write(char_log, int_fmt) info
      call log_message("LAPACK routine zgeev failed!", level="warning")
      call log_message( &
        "value for the info parameter: " // adjustl(char_log), &
        level="warning", &
        use_prefix=.false. &
      )
    end if

    deallocate(work)
    deallocate(rwork)

    call check_small_values(omega)
  end subroutine qr_invert

end submodule smod_qr_invert