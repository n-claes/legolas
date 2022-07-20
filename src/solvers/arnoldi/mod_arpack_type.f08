!> Contains a dedicated type for the various settings of the ARPACK solvers.
!! All variables that are used in different solver settings are defined, initialised
!! and set in this module.
module mod_arpack_type
  use mod_logging, only: log_message, str
  use mod_global_variables, only: dp
  implicit none

  !> General type containing the ARPACK configuration.
  type, public :: arpack_t
    !> mode of the solver
    integer, private :: mode
    !> reverse communication flag
    integer :: ido
    !> type of the matrix B ("I" = unit matrix, "G" = general matrix)
    character(len=1), private :: bmat
    !> dimension of the eigenvalue problem
    integer, private :: evpdim
    !> which eigenvalues to calculate
    character(len=2), private :: which
    !> number of eigenvalues to calculate
    integer, private :: nev
    !> stopping criteria, relative accuracy of Ritz eigenvalues
    real(dp), private :: tolerance
    !> residual vector
    complex(dp), allocatable :: residual(:)
    !> indicates how many Arnoldi vectors are generated each iteration
    integer, private :: ncv
    !> integer array containing mode and parameters
    integer :: iparam(11)
    !> length of workl array, must be at least 3 * ncv**2 + 5 * ncv
    integer, private :: lworkl
    !> info parameter
    integer :: info
    !> maximum number of iterations
    integer, private :: maxiter


    contains

    procedure, public :: get_bmat
    procedure, public :: get_evpdim
    procedure, public :: get_which
    procedure, public :: get_nev
    procedure, public :: get_tolerance
    procedure, public :: get_ncv
    procedure, public :: get_lworkl
    procedure, public :: parse_znaupd_info
    procedure, public :: parse_zneupd_info
    procedure, public :: parse_finished_stats
    procedure, public :: destroy

    procedure, private :: set_mode
    procedure, private :: set_bmat
    procedure, private :: set_which
    procedure, private :: set_nev
    procedure, private :: set_ncv
    procedure, private :: set_maxiter
  end type arpack_t

  private

  public :: new_arpack_config


contains

  !> Constructor for a new ARPACK configuration based on the dimension of the eigenvalue
  !! problem, mode of the solver and type of the B-matrix. Initialises required
  !! variables and allocates work arrays to be used when calling the solvers.
  function new_arpack_config( &
    evpdim, mode, bmat, which, nev, tolerance, maxiter, ncv &
  ) result(arpack_config)
    !> dimension of the eigenvalue problem
    integer, intent(in) :: evpdim
    !> mode for the solver
    integer, intent(in) :: mode
    !> type of the matrix B
    character(len=1), intent(in) :: bmat
    !> which eigenvalues to calculate
    character(len=2), intent(in) :: which
    !> number of eigenvalues to calculate
    integer, intent(in) :: nev
    !> relative accuracy (stopping criteria) for eigenvalues
    real(dp), intent(in) :: tolerance
    !> maximum number of iterations, defaults to 10*evpdim
    integer, intent(in) :: maxiter
    !> number of Arnoldi basis vectors
    integer, intent(in) :: ncv
    !> initialised arpack configuration
    type(arpack_t) :: arpack_config

    call log_message("configuring Arnoldi parameters", level="debug")
    arpack_config%evpdim = evpdim
    call arpack_config%set_mode(mode)

    arpack_config%ido = 0  ! 0 means first call to reverse communication interface
    call arpack_config%set_bmat(bmat)
    call arpack_config%set_which(which)
    call arpack_config%set_nev(nev)
    arpack_config%tolerance = tolerance
    allocate(arpack_config%residual(evpdim))
    call arpack_config%set_ncv(ncv)
    call arpack_config%set_maxiter(maxiter)
    ! iparam(1) = ishift = 1 means restart with shifts from Hessenberg matrix
    arpack_config%iparam(1) = 1
  end function new_arpack_config


  !> Sets the mode for the solver.
  subroutine set_mode(this, mode)
    !> type instance
    class(arpack_t), intent(inout) :: this
    !> solver mode
    integer, intent(in) :: mode
    integer, parameter :: allowed_modes(3) = [1, 2, 3]

    if (.not. any(mode == allowed_modes)) then
      call log_message( &
        "Arnoldi: mode = " // str(mode) // " is invalid, expected 1, 2 or 3", &
        level="error" &
      )
      return
    end if
    this%mode = mode
    this%iparam(7) = this%mode
    call log_message("Arnoldi: mode set to " // str(this%mode), level="debug")
  end subroutine set_mode


  !> Sets the type of B-matrix.
  subroutine set_bmat(this, bmat)
    !> type instance
    class(arpack_t), intent(inout) :: this
    !> type of B-matrix
    character, intent(in) :: bmat
    character :: allowed_bmats(2) = ["I", "G"]

    if (.not. any(bmat == allowed_bmats)) then
      call log_message( &
        "Arnoldi: bmat = " // bmat // " is invalid, expected 'I' or 'G'", &
        level="error" &
      )
      return
    end if
    this%bmat = bmat
    call log_message("Arnoldi: bmat set to " // this%bmat, level="debug")
  end subroutine set_bmat


  !> Setter for the "which" argument of ARPACK routines.
  subroutine set_which(this, which)
    !> type instance
    class(arpack_t), intent(inout) :: this
    !> which kind of eigenvalues to calculate
    character(2), intent(in) :: which
    character(2) :: allowed_which(6) = ["LM", "SM", "LR", "SR", "LI", "SI"]

    if (.not. any(which == allowed_which)) then
      call log_message( &
        "Arnoldi: which_eigenvalues = " // which &
        // " is invalid, expected one of " // str(allowed_which), &
        level="error" &
      )
      return
    end if
    this%which = which
    call log_message("Arnoldi: which set to " // this%which, level="debug")
  end subroutine set_which


  !> Setter for number of eigenvalues to calculate.
  subroutine set_nev(this, nev)
    !> type instance
    class(arpack_t), intent(inout) :: this
    !> number of eigenvalues to calculate
    integer, intent(in) :: nev

    if (nev <= 0) then
      call log_message( &
        "Arnoldi: number of eigenvalues must be >= 0 but got " // str(nev), &
        level="error" &
      )
      return
    end if
    if (nev >= this%evpdim) then
      call log_message( &
        "Arnoldi: number of eigenvalues (" // str(nev) &
        // ") >= " // "matrix size (" // str(this%evpdim) // ")", &
        level="error" &
      )
      return
    end if
    this%nev = nev
    call log_message("Arnoldi: nev set to " // str(this%nev), level="debug")
  end subroutine set_nev


  !> Setter for ncv, the number of Arnoldi basis vectors to calculate.
  !! This should satisfy 1 <= ncv - nev and ncv <= evpdim, with recommended
  !! value ncv = 2 * nev (see arpack docs).
  subroutine set_ncv(this, ncv)
    !> type instance
    class(arpack_t), intent(inout) :: this
    !> value for ncv
    integer, intent(in) :: ncv

    if (ncv == 0) then
      this%ncv = min(this%get_evpdim(), 2 * this%get_nev())
    else
      this%ncv = ncv
    end if
    if (1 > this%ncv - this%get_nev()) then
      call log_message( &
        "requesting too many eigenvalues, expected ncv - nev > 1 but got " &
        // str(this%ncv - this%get_nev()), &
        level="error" &
      )
      return
    end if
    if (this%ncv > this%get_evpdim()) then
      call log_message( &
        "ncv too high, expected ncv < N but got ncv = " // str(this%ncv) &
        // " and N = " // str(this%get_evpdim()), &
        level="error" &
      )
      return
    end if
    call log_message("Arnoldi: ncv set to " // str(this%ncv), level="debug")
  end subroutine set_ncv


  !> Sets the maximum number of iterations that ARPACK is allowed to take, defaults
  !! to 10 * N with N the dimension of the eigenvalue problem.
  !! @warning Throws a warning if <tt>maxiter</tt> is smaller than 10*N. @endwarning
  subroutine set_maxiter(this, maxiter)
    !> type instance
    class(arpack_t), intent(inout) :: this
    !> maximum number of iterations
    integer, intent(in) :: maxiter
    integer :: min_maxiter

    min_maxiter = 10 * this%get_evpdim()
    if (maxiter < 0) then
      call log_message( &
        "Arnoldi: maxiter must be positive, but is equal to " // str(maxiter), &
        level="error" &
      )
      return
    end if
    if (maxiter == 0) then
      this%maxiter = min_maxiter
    else
      this%maxiter = maxiter
    end if
    if (this%maxiter < min_maxiter) then ! LCOV_EXCL_START
      call log_message( &
        "Arnoldi: maxiter (" // str(maxiter) // ") below recommended 10*N (" &
        // str(min_maxiter) // ")", &
        level="warning" &
      )
    end if ! LCOV_EXCL_STOP
    this%iparam(3) = this%maxiter
  end subroutine set_maxiter


  !> Getter for kind of B-matrix in eigenvalue problem.
  pure character function get_bmat(this)
    !> type instance
    class(arpack_t), intent(in) :: this

    get_bmat = this%bmat
  end function get_bmat


  !> Getter for dimension of eigenvalue problem.
  pure integer function get_evpdim(this)
    !> type instance
    class(arpack_t), intent(in) :: this

    get_evpdim = this%evpdim
  end function get_evpdim


  !> Getter for which eigenvalues to return.
  pure character(2) function get_which(this)
    !> type instance
    class(arpack_t), intent(in) :: this

    get_which = this%which
  end function get_which


  !> Getter for number of eigenvalues to calculate.
  pure integer function get_nev(this)
    !> type instance
    class(arpack_t), intent(in) :: this

    get_nev = this%nev
  end function get_nev


  !> Getter for tolerance (relative accuracy) to indicate eigenvalue convergence.
  pure real(dp) function get_tolerance(this)
    !> type instance
    class(arpack_t), intent(in) :: this

    get_tolerance = this%tolerance
  end function get_tolerance


  !> Getter for number of Arnoldi basis vectors that should be calculated.
  pure integer function get_ncv(this)
    !> type instance
    class(arpack_t), intent(in) :: this

    get_ncv = this%ncv
  end function get_ncv


  !> Getter for length of workl array, returns 3 * ncv**2 + 5 * ncv
  pure integer function get_lworkl(this)
    !> type instance
    class(arpack_t), intent(in) :: this
    integer :: ncv

    ncv = this%get_ncv()
    get_lworkl = 3 * ncv**2 + 5 * ncv
  end function get_lworkl


  !> Destructor, deallocates variables.
  pure subroutine destroy(this)
    !> type instance
    class(arpack_t), intent(inout) :: this

    deallocate(this%residual)
  end subroutine destroy


  !> Parses the info parameter that comes out of ARPACK's <tt>znaupd</tt> method.
  !! If info = 0, everything behaved nicely and the reverse communication subroutines
  !! exited properly. If info is any other value something went wrong and
  !! we handle it accordingly.
  subroutine parse_znaupd_info(this, converged)
    !> reference to type object
    class(arpack_t), intent(in) :: this
    !> if .true. the reverse communication routines converged, .false. otherwise
    logical, intent(out)  :: converged

    call log_message("checking znaupd info parameter", level="debug")
    converged = .false.
    select case(this % info)
    case(0)
      converged = .true.
    case(1) ! LCOV_EXCL_START
      call log_message("ARPACK failed to converge! (maxiter reached)", level="warning")
      call log_message("number of iterations: " // str(this % maxiter), level="warning")
      call log_message( &
        "number of converged eigenvalues: " // str(this % iparam(5)) // &
        " / " // str(this % nev), &
        level="warning" &
      )
    case(3)
      call log_message( &
        "znaupd: no shifts could be applied during Arnoldi iteration", level="error" &
      )
      return
    case(-6)
      call log_message("znaupd: bmat must be 'I' or 'G'", level="error")
    case(-8)
      call log_message("znaupd: error LAPACK eigenvalue calculation", level="error")
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
    class(arpack_t), intent(in)  :: this

    call log_message("checking zneupd info parameter", level="debug")

    select case(this % info)
    case(0)
      return
    case(-8) ! LCOV_EXCL_START
      call log_message( &
        "zneupd: error from LAPACK eigenvalue calculation", level="error" &
      )
      return
    case(-9)
      call log_message( &
        "zneupd: error from LAPACK eigenvector calculation (ztrevc)", level="error" &
      )
      return
    case(-14)
      call log_message( &
        "zneupd: no eigenvalues with sufficient accuracy found", level="error" &
      )
      return
    case(-15)
      call log_message( &
        "zneupd: different count for converged eigenvalues than znaupd", level="error" &
      )
      return
    case default
      call log_message( &
        "zneupd: unexpected info = " // str(this % info) // " value", level="error" &
      )
      return ! LCOV_EXCL_STOP
    end select
  end subroutine parse_zneupd_info


  !> Parses the statistics that come out of ARPACK when the run is finished. Displays
  !! the number of OP*X and B*X operations and the number of re-orthogonalisation
  !! steps that were needed.
  subroutine parse_finished_stats(this)
    !> type instance
    class(arpack_t), intent(in) :: this

    call log_message("Arnoldi iteration finished. Statistics: ", level="info")
    call log_message( &
      "   Total number of OP*x operations: " // str(this%iparam(9)), &
      level="info", &
      use_prefix=.false. &
    )
    call log_message( &
      "   Total number of B*x operations: " // str(this%iparam(10)), &
      level="info", &
      use_prefix=.false. &
    )
    call log_message( &
      "   Total number of re-orthogonalisation steps: " // str(this%iparam(11)), &
      level="info", &
      use_prefix=.false. &
    )
  end subroutine parse_finished_stats
end module mod_arpack_type
