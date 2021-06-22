! =============================================================================
!> Module containing a dedicated type for the ARPACK solvers.
!! All needed variables are defined and allocated here, along with sanity
!! checks when setting certain values.
module mod_arpack_type
  use mod_logging, only: log_message, str
  use mod_global_variables, only: dp, dp_LIMIT

  implicit none

  type arpack_type
    !> reverse communication flag
    integer           :: ido
    !> maximum number of iterations
    integer           :: maxiter
    !> mode for the solver
    integer           :: mode
    !> specifies type of B-matrix ("I" = unit matrix, "G" = general)
    character(len=1)  :: bmat
    !> dimension of the eigenvalue problem
    integer           :: evpdim
    !> which eigenvalues to calculate
    character(len=2)  :: which
    !> number of eigenvalues to calculate
    integer           :: nev
    !> stopping criteria, relative accuracy of Ritz eigenvalues
    real(dp)          :: tol
    !> residual vector used for initialisation
    complex(dp), allocatable  :: residual(:)
    !> number of Arnoldi basis vectors
    integer                   :: ncv
    !> contains the Arnoldi basis vectors
    complex(dp), allocatable  :: arnoldi_vectors(:, :)
    !> integer array containing mode and parameters
    integer                   :: iparam(11)
    !> integer array containing pointers to mark work array locations
    integer                   :: ipntr(14)
    !> complex work array, length 3N
    complex(dp), allocatable  :: workd(:)
    !> complex work array, length lworkl
    complex(dp), allocatable  :: workl(:)
    !> length of workl array, at least 3*ncv**2 + 5*ncv
    integer                   :: lworkl
    !> real work array, length ncv
    real(dp), allocatable     :: rwork(:)
    !> info parameter
    integer                   :: info
    !> if .true., also calculate eigenvectors
    logical                   :: rvec
    !> specifies form of basis: "A" = nev Ritz vectors, "P" = nev Schur vectors
    character(len=1)          :: howmny
    !> logical array of dimension ncv, selects Ritz vectors
    logical, allocatable      :: select_vectors(:)
    !> represents the shift (only referenced if shift-invert is used)
    complex(dp)               :: sigma
    !> work array for eigenvalues, size 2*ncv
    complex(dp), allocatable  :: workev(:)
    !> number of converged eigenvalues
    integer                   :: nconv

    contains

      !> initialises derived type
      procedure, public   :: initialise
      !> sets mode of ARPACK solver
      procedure, public   :: set_mode
      !> sets sigma for shift-invert mode
      procedure, public   :: set_sigma
      !> parses info parameter of <tt>znaupd</tt>
      procedure, public   :: parse_znaupd_info
      !> parses info parameter of <tt>zneupd</tt>
      procedure, public   :: parse_zneupd_info
      !> cleanup routine
      procedure, public   :: tear_down
      !> setter for requested number of eigenvalues
      procedure, private  :: set_nev
      !> setter for 'which' argument
      procedure, private  :: set_which
      !> setter for maximum number of iterations to take
      procedure, private  :: set_maxiter
  end type arpack_type

  private

  public :: arpack_type

contains

  !> Initialises the derived type based on the given dimension
  !! of the eigenvalue problem. Everything is allocated and the variables
  !! (<tt>nev, which, maxiter, iparam, ido</tt> etc.)
  !! are prepared for their calls in the ARPACK routines.
  subroutine initialise(this, evpdim)
    !> reference to type object
    class(arpack_type)  :: this
    !> dimension of eigenvalue problem
    integer, intent(in) :: evpdim

    ! set dimension and related parameters
    this % evpdim = evpdim
    call this % set_nev()
    call this % set_which()
    call this % set_maxiter()
    this % ncv = min(this % evpdim, 2 * this % nev)
    this % tol = dp_LIMIT

    ! allocate and initialise Arnoldi vectors
    allocate(this % arnoldi_vectors(this % evpdim, this % ncv))
    this % arnoldi_vectors = (0.0d0, 0.0d0)

    ! allocate work arrays
    this % lworkl = 3 * this % ncv * (this % ncv + 2)
    allocate(this % workl(this % lworkl))
    allocate(this % rwork(this % ncv))
    allocate(this % workd(3 * this % evpdim))

    ! set parameters and solver mode
    this % iparam(1) = 1    ! select implicit shifts (1 = restart)
    this % iparam(3) = this % maxiter   ! maximum number of iterations
    this % iparam(4) = 1    ! blocksize, HAS to be 1
    this % ido = 0    ! 0 means first call to interface
    this % info = 0   ! info = 0 at input means use a random residual starting vector
    allocate(this % residual(this % evpdim))

    ! allocate eigenvalue-related things
    allocate(this % select_vectors(this % ncv))
    allocate(this % workev(2 * this % ncv))
    this % rvec = .true.    ! always calculate eigenvectors, not expensive in ARPACK
    this % howmny = "A"     ! currently hardcoded to Ritz vectors
  end subroutine initialise


  !> Sets the mode for the ARPACK solver and passes this on
  !! to the <tt>iparam</tt> array.
  !! @warning Throws an error if mode is not equal to 1, 2 or 3. @endwarning
  subroutine set_mode(this, mode)
    !> reference to type object
    class(arpack_type)  :: this
    !> the mode to set, should be 1, 2 or 3
    integer, intent(in) :: mode

    if (mode < 1 .or. mode > 3) then
      call log_message( &
        "mode must be 1, 2 or 3 but mode = " // str(mode) // " was given", &
        level="error" &
      )
      return
    end if
    this % mode = mode
    this % iparam(7) = this % mode
  end subroutine set_mode


  !> Sets the sigma value for the shift-invert mode of ARPACK.
  !! Sigma can't be zero since the A-matrix can have zero rows, and then
  !! we run into troubles.
  !! @warning Throws an error if sigma = 0. @endwarning
  subroutine set_sigma(this, sigma)
    use mod_check_values, only: is_equal

    !> reference to type object
    class(arpack_type)      :: this
    !> sigma for the shift-invert method
    complex(dp), intent(in) :: sigma

    if (is_equal(sigma, (0.0d0, 0.0d0))) then
      call log_message( &
        "ARPACK shift-invert: sigma can not be equal to zero", &
        level="error" &
      )
      return
    end if
    this % sigma = sigma
  end subroutine set_sigma


  !> Parses the info parameter that comes out of ARPACK's <tt>znaupd</tt> method.
  !! If info = 0, everything behaved nicely the reverse communication subroutines
  !! exited properly. If info is any other value something went wrong and
  !! we handle it accordingly.
  subroutine parse_znaupd_info(this, converged)
    !> reference to type object
    class(arpack_type)    :: this
    !> if .true. the reverse communication routines converged, .false. otherwise
    logical, intent(out)  :: converged

    converged = .false.

    select case(this % info)
    case(0)
      converged = .true.
    case(1) ! LCOV_EXCL_START
      call log_message("ARPACK failed to converge! (maxiter reached)", level="warning")
      call log_message( &
        "number of iterations: " // str(this % maxiter), &
        level="warning", &
        use_prefix=.false. &
      )
      call log_message( &
        "number of converged eigenvalues: " // str(this % iparam(5)) // &
        " / " // str(this % nev), &
        level="warning", &
        use_prefix=.false. &
      )
    case(3)
      call log_message( &
        "znaupd: no shifts could be applied during Arnoldi iteration", &
        level="error" &
      )
      return
    case(-6)
      call log_message("znaupd: bmat must be 'I' or 'G'", level="error")
    case(-8)
      call log_message( &
        "znaupd: error from LAPACK eigenvalue calculation", &
        level="error" &
      )
      return
    case(-11)
      call log_message("mode = 1 and bmat = 'G' are incompatible", level="error")
      return
    case(-9999)
      call log_message("ARPACK could not build, something went wrong", level="error")
      return
    case default
      call log_message( &
        "znaupd: unexpected info = " // str(this % info) // " encountered", &
        level="error" &
      )
      return ! LCOV_EXCL_STOP
    end select
  end subroutine parse_znaupd_info


  !> Parses the info parameter that comes out of ARPACK's <tt>zneupd</tt> method.
  !! If info = 0, the eigenvalues extraction routines exited properly, if info
  !! is any other value something went wrong and we handle it accordingly.
  subroutine parse_zneupd_info(this)
    !> reference to type object
    class(arpack_type)  :: this

    select case(this % info)
    case(0)
      return
    case(-8) ! LCOV_EXCL_START
      call log_message( &
        "zneupd: error from LAPACK eigenvalue calculation", &
        level="error" &
      )
      return
    case(-9)
      call log_message( &
        "zneupd: error from LAPACK eigenvector calculation (ztrevc)", &
        level="error" &
      )
      return
    case(-14)
      call log_message( &
        "zneupd: no eigenvalues with sufficient accuracy found", &
        level="error" &
      )
      return
    case(-15)
      call log_message( &
        "zneupd: different count for converged eigenvalues than znaupd", &
        level="error" &
      )
      return
    case default
      call log_message( &
        "zneupd: unexpected info = " // str(this % info) // " value", &
        level="error" &
      )
      return ! LCOV_EXCL_STOP
    end select
  end subroutine parse_zneupd_info


  !> Setter for the number of eigenvalues that should be calculated.
  !! The requested number of eigenvalues should be positive and smaller than
  !! the dimension of the eigenvalue problem.
  subroutine set_nev(this)
    use mod_global_variables, only: number_of_eigenvalues

    !> reference to type object
    class(arpack_type)  :: this

    if (number_of_eigenvalues <= 0) then
      call log_message( &
        "number_of_eigenvalues must be >= 0, but is equal to " &
        // str(number_of_eigenvalues), &
        level="error" &
      )
      return
    end if
    if (number_of_eigenvalues >= this % evpdim) then
      call log_message( &
        "number_of_eigenvalues larger than matrix size! (" &
        // str(number_of_eigenvalues) // " > " // str(this % evpdim) // ")", &
        level="error" &
      )
      return
    end if
    this % nev = number_of_eigenvalues
  end subroutine set_nev


  !> Setter for the 'which' argument of ARPACK routines.
  !! Should be one of the following: "LM", "SM", "LR", "SR", "LI" or "SI".
  subroutine set_which(this)
    use mod_global_variables, only: which_eigenvalues

    class(arpack_type)  :: this
    character(2)        :: allowed_which(6) = ["LM", "SM", "LR", "SR", "LI", "SI"]

    if (.not. any(which_eigenvalues == allowed_which)) then
      call log_message( &
        "which_eigenvalues = " // which_eigenvalues // " is invalid", &
        level="error" &
      )
      return
    end if
    this % which = which_eigenvalues
  end subroutine set_which


  !> Setter for the maximum number of iterations that ARPACK can take,
  !! defaults to 10N with N the dimension of the eigenvalue problem.
  !! @warning Throws a warning if <tt>maxiter</tt> is smaller than 10*N. @endwarning
  subroutine set_maxiter(this)
    use mod_global_variables, only: maxiter

    class(arpack_type)  :: this

    ! is maxiter is not set in the parfile it's still 0, default to 10*N
    if (maxiter == 0) then
      maxiter = 10 * this % evpdim
    else if (maxiter < 0) then
      call log_message( &
        "maxiter has to be positive, but is equal to " // str(maxiter), level="error" &
      )
      return
    else if (maxiter < 10 * this % evpdim) then ! LCOV_EXCL_START
      call log_message( &
        "maxiter is below recommended 10*N: (" &
        // str(maxiter) // " < " // str(10 * this % evpdim) // ")", &
        level="warning" &
      )
    end if ! LCOV_EXCL_STOP
    this % maxiter = maxiter
  end subroutine set_maxiter


  !> Cleanup routine, deallocates type attributes.
  subroutine tear_down(this)
    !> reference to type object
    class(arpack_type)  :: this

    if (allocated(this % arnoldi_vectors)) then
      deallocate(this % arnoldi_vectors)
    end if
    if (allocated(this % workl)) then
      deallocate(this % workl)
    end if
    if (allocated(this % rwork)) then
      deallocate(this % rwork)
    end if
    if (allocated(this % workd)) then
      deallocate(this % workd)
    end if
    if (allocated(this % residual)) then
      deallocate(this % residual)
    end if
    if (allocated(this % select_vectors)) then
      deallocate(this % select_vectors)
    end if
    if (allocated(this % workev)) then
      deallocate(this % workev)
    end if
  end subroutine tear_down

end module mod_arpack_type
