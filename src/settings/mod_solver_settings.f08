module mod_solver_settings
  use mod_global_variables, only: dp
  implicit none

  private

  type, public :: solvers_t
    character(:), allocatable, private :: solver
    character(:), allocatable, private :: arpack_mode
    integer :: number_of_eigenvalues
    character(len=2) :: which_eigenvalues
    integer :: maxiter
    complex(dp) :: sigma
    integer :: ncv
    real(dp) :: tolerance

  contains

    procedure, public :: set_solver
    procedure, public :: get_solver
    procedure, public :: set_arpack_mode
    procedure, public :: get_arpack_mode
    procedure, public :: delete
  end type solvers_t

  public :: new_solver_settings

contains

  pure function new_solver_settings() result(solvers)
    use mod_global_variables, only: dp_LIMIT, NaN

    type(solvers_t) :: solvers

    call solvers%set_solver("QR-invert")
    call solvers%set_arpack_mode("standard")
    solvers%number_of_eigenvalues = 10
    solvers%which_eigenvalues = "LM"
    ! these two get determined at runtime
    solvers%maxiter = 0
    solvers%ncv = 0
    solvers%sigma = cmplx(NaN, NaN, kind=dp)
    solvers%tolerance = dp_LIMIT
  end function new_solver_settings


  pure subroutine set_solver(this, solver)
    class(solvers_t), intent(inout) :: this
    character(len=*), intent(in) :: solver

    this%solver = trim(adjustl(solver))
  end subroutine set_solver


  pure function get_solver(this) result(solver)
    class(solvers_t), intent(in) :: this
    character(:), allocatable :: solver

    solver = trim(adjustl(this%solver))
  end function get_solver


  pure subroutine set_arpack_mode(this, arpack_mode)
    class(solvers_t), intent(inout) :: this
    character(len=*), intent(in) :: arpack_mode

    this%arpack_mode = trim(adjustl(arpack_mode))
  end subroutine set_arpack_mode


  pure function get_arpack_mode(this) result(arpack_mode)
    class(solvers_t), intent(in) :: this
    character(:), allocatable :: arpack_mode

    arpack_mode = trim(adjustl(this%arpack_mode))
  end function get_arpack_mode


  pure subroutine delete(this)
    class(solvers_t), intent(inout) :: this

    if (allocated(this%solver)) deallocate(this%solver)
    if (allocated(this%arpack_mode)) deallocate(this%arpack_mode)
  end subroutine delete

end module mod_solver_settings
