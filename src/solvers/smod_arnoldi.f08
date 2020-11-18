submodule (mod_solvers) smod_arnoldi
  use mod_matrix_operations, only: invert_matrix, multiply_matrices
  use mod_arpack_type, only: arpack_type
  use mod_global_variables, only: arpack_mode
  implicit none

contains


  module subroutine arnoldi(matrix_A, matrix_B, omega, vr)
    !> matrix A
    complex(dp), intent(in)   :: matrix_A(:, :)
    !> matrix B
    real(dp), intent(in)      :: matrix_B(:, :)
    !> array with calculated eigenvalues
    complex(dp), intent(out)  :: omega(:)
    !> array with right eigenvectors
    complex(dp), intent(out)  :: vr(:, :)

#if _ARPACK_FOUND
    !> inverse B-matrix
    real(dp), allocatable    :: B_inv(:, :)
    !> matrix \(B^{-1}A\)
    complex(dp), allocatable :: B_invA(:, :)
    !> matrix product A*x
    complex(dp) :: ax(size(matrix_A, dim=1))
    !> dedicated type with all ARPACK-related stuff
    type(arpack_type) :: arpackparams
    !> flag to check if solver converged
    logical, save     :: converged

    integer   :: i, xleft, xright, yleft, yright

    ! initialise arpack params
    call arpackparams % initialise(evpdim=size(matrix_A, dim=1))

    ! cycle through possible modes
    select case(arpack_mode)
    case("standard")
      ! solves standard eigenvalue problem Ax = wx
      allocate(B_inv(size(matrix_B, dim=1), size(matrix_B, dim=2)))
      allocate(B_invA(size(matrix_A, dim=1), size(matrix_A, dim=2)))
      ! do inversion of B
      call invert_matrix(matrix_B, B_inv)
      ! do matrix multiplication B^{-1}A
      call multiply_matrices(B_inv, matrix_A, B_invA)
      deallocate(B_inv)   ! no longer used after this

      call arpackparams % set_mode(1)
      ! in this case B is the identity matrix
      arpackparams % bmat = "I"
    case default
      call arpackparams % tear_down()
      call log_message("unknown mode for ARPACK: " // arpack_mode, level="error")
      return
    end select

    converged = .false.
    ! keep iterating as long as the eigenvalues are not converged.
    ! if convergence is achieved or the maximum number of iterations is reached,
    ! ARPACK sets ido=99 so while loop breaks, we know what happened through 'info'.
    ! TODO: if mode = 2 or mode = 3 we can also have ido = 2 and ido = 3
    do while (.not. converged)
      call znaupd( &
        arpackparams % ido, &
        arpackparams % bmat, &
        arpackparams % evpdim, &
        arpackparams % which, &
        arpackparams % nev, &
        arpackparams % tol, &
        arpackparams % residual, &
        arpackparams % ncv, &
        arpackparams % arnoldi_vectors, &
        size(arpackparams % arnoldi_vectors, dim=1), &
        arpackparams % iparam, &
        arpackparams % ipntr, &
        arpackparams % workd, &
        arpackparams % workl, &
        arpackparams % lworkl, &
        arpackparams % rwork, &
        arpackparams % info &
      )
      if (arpackparams % ido == -1 .or. arpackparams % ido == 1) then
        ! do matrix vector multiplication A*x -> y
        ! start of x is given by workd(ipntr(1))
        xleft = arpackparams % ipntr(1)
        xright = xleft + arpackparams % evpdim - 1
        yleft = arpackparams % ipntr(2)
        yright = yleft + arpackparams % evpdim - 1
        ! start of y is stored in workd(ipntr(2))
        call multiply_matrices( &
          B_invA, &
          arpackparams % workd(xleft:xright), &
          arpackparams % workd(yleft:yright) &
        )
      else
        exit
      end if
    end do

    ! check info parameter from znaupd, this errors if necessary
    call arpackparams % parse_znaupd_info(converged)

    ! if we have a normal exit, extract the eigenvalues through zneupd
    call zneupd( &
      arpackparams % rvec, &
      arpackparams % howmny, &
      arpackparams % select_vectors, &
      omega(1:arpackparams % nev), &
      vr(:, 1:arpackparams % nev), &
      size(vr, dim=1), &
      arpackparams % sigma, &
      arpackparams % workev, &
      arpackparams % bmat, &
      arpackparams % evpdim, &
      arpackparams % which, &
      arpackparams % nev, &
      arpackparams % tol, &
      arpackparams % residual, &
      arpackparams % ncv, &
      arpackparams % arnoldi_vectors, &
      size(arpackparams % arnoldi_vectors, dim=1), &
      arpackparams % iparam, &
      arpackparams % ipntr, &
      arpackparams % workd, &
      arpackparams % workl, &
      arpackparams % lworkl, &
      arpackparams % rwork, &
      arpackparams % info &
    )

    ! check info parameter from zneupd, this errors if necessary
    call arpackparams % parse_zneupd_info()

    ! calculate residual vector || A*x - lambda*x ||
    arpackparams % nconv = arpackparams % iparam(5)
    allocate(residual_norm(arpackparams % nconv))    ! defined in parent module
    do i = 1, arpackparams % nconv
      call multiply_matrices(B_invA, vr(:, i), ax)
      residual_norm(i) = real(sqrt(sum( &
        (ax - omega(i) * vr(:, i)) * conjg(ax - omega(i) * vr(:, i)) &
      )))
    end do

    call arpackparams % tear_down()

#else
    call log_message( &
      "ARPACK was not found and/or CMake failed to link", &
      level="warning" &
    )
    call log_message( &
      "unable to use 'arnoldi, try another solver!", &
      level="error" &
    )
#endif

  end subroutine arnoldi

end submodule smod_arnoldi