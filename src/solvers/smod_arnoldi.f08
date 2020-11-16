submodule (mod_solvers) smod_arnoldi
  use mod_global_variables, only: arpack_mode, number_of_eigenvalues, &
    which_eigenvalues, maxiter, dp_LIMIT
  use mod_matrix_operations, only: invert_matrix, multiply_matrices
  use mod_logging, only: char_log2
  implicit none

  complex(dp), allocatable  :: residual(:)

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
    ! make these allocatable, so we only use memory if we are actually using them
    !> inverse B-matrix
    real(dp), allocatable    :: B_inv(:, :)
    !> matrix \(B^{-1}A\)
    complex(dp), allocatable :: B_invA(:, :)

    integer            :: mode
    !> type of the matrix B, "I" means unit matrix, "G" means general matrix
    character(len=1)   :: bmat
    !> dimension of the eigenvalue problem
    integer            :: evpdim
    !> number of Arnoldi vectors generated
    integer            :: ncv
    !> contains the Arnoldi basis vectors
    complex(dp), allocatable  :: arnoldi_vectors(:, :)
    !> integer array containing mode and parameters
    integer            :: iparam(11)
    !> integer array containing pointers to mark work array locations
    integer            :: ipntr(14)
    !> work array of length 3*N
    complex(dp)        :: workd(3 * size(matrix_A, dim=1))
    !> length of lworkl array
    integer            :: lworkl
    !> lworkl array, needs size of at least 3*ncv**2 + 5*ncv
    complex(dp), allocatable  :: workl(:)
    !> rwork array of length ncv
    real(dp), allocatable     :: rwork(:)
    !> info parameter
    integer            :: info
    !> reverse communication flag
    integer            :: ido
    !> flag to check if solver converged
    logical, save      :: converged
    !> logical array to select Ritz vectors, size ncv
    logical, allocatable :: select_vectors(:)
    !> represents the shift
    complex(dp)        :: sigma
    !> work array for eigenvalues, size 2*ncv
    complex(dp), allocatable  :: workev(:)

    evpdim = size(matrix_A, dim=1)
    allocate(residual(evpdim))

    ! is maxiter is not set in the parfile it's still 0, default to 10*N
    if (maxiter == 0) then
      maxiter = 10 * evpdim
    end if
    residual = (0.0d0, 0.0d0)
    ! TODO: allow for customisation of ncv through parfile
    ncv = min(evpdim, 2 * number_of_eigenvalues)
    allocate(arnoldi_vectors(evpdim, ncv))
    arnoldi_vectors = (0.0d0, 0.0d0)
    lworkl = 3 * ncv * (ncv + 2)
    allocate(workl(lworkl))
    allocate(rwork(ncv))

    call do_arpack_sanity_checks(evpdim=evpdim)

    ! cycle through possible modes
    select case(arpack_mode)
    case("standard")
      ! solves standard eigenvalue problem Ax = wx
      mode = 1

      allocate(B_inv(size(matrix_B, dim=1), size(matrix_B, dim=2)))
      allocate(B_invA(size(matrix_A, dim=1), size(matrix_A, dim=2)))
      ! do inversion of B
      call invert_matrix(matrix_B, B_inv)
      ! do matrix multiplication B^{-1}A
      call multiply_matrices(B_inv, matrix_A, B_invA)
      deallocate(B_inv)   ! no longer used after this

      bmat = "I"
    case default
      call log_message("unknown mode for ARPACK: " // arpack_mode, level="error")
      return
    end select

    iparam(1) = 1           ! select implicit shifts (1 = restart)
    iparam(3) = maxiter     ! maximum number of iterations
    iparam(4) = 1           ! blocksize, HAS to be 1
    iparam(7) = mode        ! mode for the solver

    ido = 0                 ! 0 means first call to interface
    converged = .false.

    ! keep iterating as long as the eigenvalues are not converged.
    ! if convergence is achieved or the maximum number of iterations is reached,
    ! ARPACK sets ido=99 so while loop breaks, we know what happened through 'info'.
    ! TODO: if mode = 2 or mode = 3 we can also have ido = 2 and ido = 3
    do while (.not. converged)
      call znaupd( &
        ido, &
        bmat, &
        evpdim, &
        which_eigenvalues, &
        number_of_eigenvalues, &
        dp_LIMIT, &
        residual, &
        ncv, &
        arnoldi_vectors, &
        size(arnoldi_vectors, dim=1), &
        iparam, &
        ipntr, &
        workd, &
        workl, &
        lworkl, &
        rwork, &
        info &
      )
      if (ido == -1 .or. ido == 1) then
        ! do matrix vector multiplication A*x -> y
        ! start of x is given by workd(ipntr(1))
        ! start of y is stored in workd(ipntr(2))
        call multiply_matrices( &
          B_invA, &
          workd(ipntr(1):ipntr(1) + evpdim - 1), &
          workd(ipntr(2):ipntr(2) + evpdim - 1) &
        )
      else
        exit
      end if
    end do

    ! find out what happened on exit
    if (info == 0) then
      ! info = 0 means normal exit
      converged = .true.
    else if (info == 1) then
      ! info = 1 means maximum number of iterations reached
      call log_message( &
        "ARPACK failed to converge! (maxiter reached)", &
        level="warning" &
      )
      write(char_log, int_fmt) maxiter
      call log_message( &
        "number of iterations: " // trim(adjustl(char_log)), &
        level="warning", &
        use_prefix=.false. &
      )
      write(char_log, int_fmt) iparam(5)
      write(char_log2, int_fmt) number_of_eigenvalues
      call log_message( &
        "number of converged eigenvalues: " // trim(adjustl(char_log)) // &
          " / " // trim(adjustl(char_log2)), &
        level="warning", &
        use_prefix=.false. &
      )
    else
      write(char_log, int_fmt) info
      call log_message( &
        "znaupd: unexpected info = " // trim(adjustl(char_log)) // " value", &
        level="error" &
      )
      return
    end if

    ! if we have a normal exit, extract the eigenvalues through zneupd
    omega = (0.0d0, 0.0d0)
    allocate(select_vectors(ncv))
    allocate(workev(2 * ncv))
    call zneupd( &
      write_eigenfunctions, &
      "A", &
      select_vectors, &
      omega(1:number_of_eigenvalues + 1), &
      vr(:, 1:number_of_eigenvalues), &
      size(vr, dim=1), &
      sigma, &
      workev, &
      bmat, &
      evpdim, &
      which_eigenvalues, &
      number_of_eigenvalues, &
      dp_LIMIT, &
      residual, &
      ncv, &
      arnoldi_vectors, &
      size(arnoldi_vectors, dim=1), &
      iparam, &
      ipntr, &
      workd, &
      workl, &
      lworkl, &
      rwork, &
      info &
    )

    if (info /= 0) then
      write(char_log, int_fmt) info
      call log_message( &
        "zneupd: unexpected info = " // trim(adjustl(char_log)) // " value", &
        level="error" &
      )
      return
    end if

    if (allocated(B_invA)) then
      deallocate(B_invA)
    end if
    deallocate(arnoldi_vectors)
    deallocate(workl)
    deallocate(rwork)
    deallocate(select_vectors)
    deallocate(workev)

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


  subroutine do_arpack_sanity_checks(evpdim)
    integer, intent(in) :: evpdim
    character(2)        :: allowed_which(6)

    ! validity check for requested number of eigenvalues
    if (number_of_eigenvalues <= 0) then
      write(char_log, int_fmt) number_of_eigenvalues
      call log_message( &
        "number_of_eigenvalues must be >= 0, but is equal to " &
          // adjustl(trim(char_log)), &
        level="error" &
      )
      return
    end if
    if (number_of_eigenvalues >= evpdim) then
      write(char_log, int_fmt) number_of_eigenvalues
      write(char_log2, int_fmt) evpdim
      call log_message( &
        "number_of_eigenvalues larger than matrix size! (" &
          // trim(adjustl(char_log)) // " > " // trim(adjustl(char_log2)) // ")", &
        level="error" &
      )
      return
    end if

    ! validity check for max Arnoldi iterations
    if (maxiter <= 0) then
      write(char_log, int_fmt) maxiter
      call log_message( &
        "maxiter has to be positive, but is equal to " &
          // trim(adjustl(char_log)), level="error" &
      )
      return
    else if (maxiter < 10 * evpdim) then
      write(char_log, int_fmt) maxiter
      write(char_log2, int_fmt) 10 * evpdim
      call log_message( &
        "maxiter is below recommended 10*N: (" &
          // trim(adjustl(char_log)) // " < " // trim(adjustl(char_log2)) // ")", &
        level="warning" &
      )
    end if

    ! validity check for provided 'which' argument
    allowed_which = [character(len=2) :: "LM", "SM", "LR", "SR", "LI", "SI"]
    if (.not. any(which_eigenvalues == allowed_which)) then
      call log_message( &
        "which_eigenvalues = " // which_eigenvalues // " is invalid; must be one of &
        &the following: 'LM', 'SM', 'LR', 'SR', 'LI' or 'SI'", level="error" &
      )
      return
    end if

  end subroutine do_arpack_sanity_checks

end submodule smod_arnoldi