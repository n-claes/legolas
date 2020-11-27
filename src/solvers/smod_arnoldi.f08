! =============================================================================
!> Submodule containing the implementation of the various ARPACK solvers.
!! Preprocessor directives are used here, if ARPACK was not found during compilation
!! most of this submodule will not be compiled to avoid unknown calls to methods.
!! We first cycle through the possible modes, set up the relevant matrix operators and
!! call the corresponding methods. Then the reverse communication interface
!! <tt>znaupd</tt> is called, which will keep iterating until either
!! everything is converged or the maximum number of iterations <tt>maxiter</tt>
!! is reached. Then eigenvalues are extracted using <tt>zneupd</tt>, eigenvectors
!! are calculated and the residual of each eigenvalue is computed
!! which should be small for converged values.
submodule (mod_solvers) smod_arnoldi
  use mod_matrix_operations, only: invert_matrix, multiply_matrices
  use mod_arpack_type, only: arpack_type
  use mod_global_variables, only: arpack_mode, sigma
  implicit none

contains


  !> Implementation of the ARPACK solvers
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
    !> operator to use in the various cases
    complex(dp), allocatable :: OP(:, :)
    !> matrix product A*x used during residual calculation
    complex(dp), allocatable :: ax(:)
    !> matrix product B*x used during residual calculation
    complex(dp), allocatable :: bx(:)
    !> dedicated type with all ARPACK-related stuff
    type(arpack_type) :: arpackparams
    !> flag to check if solver converged
    logical, save     :: converged
    !> flag if we are in shift-invert mode or not
    logical, save     :: shift_invert

    integer   :: i, xleft, xright, yleft, yright

    ! initialise arpack params
    call arpackparams % initialise(evpdim=size(matrix_A, dim=1))
    shift_invert = .false.
    allocate(OP, mold=matrix_A)

    ! cycle through possible modes
    select case(arpack_mode)
    case("standard")
      call log_message("initialising Arnoldi iteration, standard mode", level="debug")
      ! solves standard eigenvalue problem Ax = wx
      call arpackparams % set_mode(1)
      ! in this case B is the identity matrix
      arpackparams % bmat = "I"
    case("general")
      call log_message("initialising Arnoldi iteration, general mode", level="debug")
      ! solves general eigenvalue problem Ax = wBx
      call arpackparams % set_mode(2)
      ! in this case B is a general matrix
      arpackparams % bmat = "G"
    case("shift-invert")
      call log_message( &
        "initialising Arnoldi iteration, shift-invert mode", &
        level="debug" &
      )
      ! solves eigenvalue problem Ax = wBx in shift-invert mode
      call arpackparams % set_mode(3)
      arpackparams % bmat = "G"
      call arpackparams % set_sigma(sigma)
      shift_invert = .true.
    case default
      call arpackparams % tear_down()
      call log_message("unknown mode for ARPACK: " // arpack_mode, level="error")
      return
    end select

    if (shift_invert) then
      ! calculate operator inv[A - sigma * B]
      call invert_matrix(matrix_A - sigma * matrix_B, OP)
    else
      ! get inverse of B matrix, needed in both standard and general eigenvalue problem
      allocate(B_inv, mold=matrix_B)
      ! do inversion of B
      call invert_matrix(matrix_B, B_inv)
      ! set operator OP <- inv[B]*A
      call multiply_matrices(B_inv, matrix_A, OP)
      deallocate(B_inv)   ! no longer used after this
    end if

    call log_message("doing Arnoldi iteration...", level="debug")
    converged = .false.
    ! keep iterating as long as the eigenvalues are not converged.
    ! if convergence is achieved or the maximum number of iterations is reached,
    ! ARPACK sets ido=99 so while loop breaks, we know what happened through 'info'.
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

      ! start of x is given by ipntr(1)
      xleft = arpackparams % ipntr(1)
      xright = xleft + arpackparams % evpdim - 1
      ! start of y is given by ipntr(2)
      yleft = arpackparams % ipntr(2)
      yright = yleft + arpackparams % evpdim - 1

      ! take action depending on "ido" outcome
      select case(arpackparams % ido)
        ! x-values are given by workd(xleft:xright)
        ! y-values are given by workd(yleft:yright)
      case(-1)
        ! entered on initialisation, forces starting vector in OP range
        call multiply_matrices( &
          OP, &
          arpackparams % workd(xleft:xright), &
          arpackparams % workd(yleft:yright) &
        )
      case(1)
        ! entered during iteration, calculates OP * x --> y
        ! mode = 1: OP = A          (but in our case A = inv[B] * A)
        ! mode = 2: OP = inv[B] * A
        ! mode = 3: OP = inv[A - sigma * B]
        if (shift_invert) then
          ! in shift-invert mode start of B*x has been saved in ipntr(3)
          ! this takes [...] * B into account
          xleft = arpackparams % ipntr(3)
          xright = xleft + arpackparams % evpdim - 1
        end if
        call multiply_matrices( &
          OP, &
          arpackparams % workd(xleft:xright), &
          arpackparams % workd(yleft:yright) &
        )
      case(2)
        ! entered during iteration, calculates B * x --> y
        call multiply_matrices( &
          matrix_B, &
          arpackparams % workd(xleft:xright), &
          arpackparams % workd(yleft:yright) &
        )
      case default
        ! entered if values converged or maxiter is reached: exit while loop
        exit
      end select
    end do

    ! check info parameter from znaupd, this errors if necessary
    call log_message("checking znaupd info parameter", level="debug")
    call arpackparams % parse_znaupd_info(converged)

    call log_message("extracting eigenvalues...", level="debug")
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
    call log_message("checking zneupd info parameter", level="debug")
    call arpackparams % parse_zneupd_info()

    ! calculate residual vector
    ! - for standard problem: || A*x - lambda * x ||
    ! - for general problem : || A*x - lambda * B * x ||
    arpackparams % nconv = arpackparams % iparam(5)
    allocate(residual_norm(arpackparams % nconv))    ! defined in parent module
    allocate(ax(arpackparams % evpdim))
    allocate(bx(arpackparams % evpdim))
    call log_message("calculating residual norm", level="debug")
    do i = 1, arpackparams % nconv
      call multiply_matrices(matrix_A, vr(:, i), ax)
      if (arpackparams % bmat == "G") then
        call multiply_matrices(matrix_B, vr(:, i), bx)
      else
        bx = vr(:, i)
      end if
      residual_norm(i) = real(sqrt(sum( &
        (ax - omega(i) * bx) * conjg(ax - omega(i) * bx) &
      )))
    end do

    deallocate(OP)
    deallocate(ax)
    deallocate(bx)
    call arpackparams % tear_down()

#else
    call log_message( &
      "ARPACK was not found and/or CMake failed to link", &
      level="warning" &
    )
    call log_message( &
      "unable to use 'arnoldi', try another solver!", &
      level="error" &
    )
#endif

  end subroutine arnoldi

end submodule smod_arnoldi