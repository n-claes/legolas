! =============================================================================
!> Parent module for everything solver-related. Interfaces to the different
!! submodules are defined here, and the <tt>solve_evp</tt> routine calls the
!! correct solver based on parfile settings.
module mod_solvers
  use mod_global_variables, only: dp, write_eigenfunctions
  use mod_logging, only: log_message, str
  use mod_check_values, only: set_small_values_to_zero
  implicit none

  !> residual norm || Ax - \(\lambda\)Bx || of the eigenvalue problem
  real(dp), allocatable  :: residual_norm(:)

  private

  !> interface to the different solution methods implemented in submodules
  interface
    module subroutine qr_invert(matrix_A, matrix_B, omega, vr)
      !> matrix A
      complex(dp), intent(in)   :: matrix_A(:, :)
      !> matrix B
      real(dp), intent(in)      :: matrix_B(:, :)
      !> array with eigenvalues
      complex(dp), intent(out)  :: omega(:)
      !> array with right eigenvectors
      complex(dp), intent(out)  :: vr(:, :)
    end subroutine

    module subroutine qz_direct(matrix_A, matrix_B, omega, vr)
      !> matrix A
      complex(dp), intent(in)   :: matrix_A(:, :)
      !> matrix B
      real(dp), intent(in)      :: matrix_B(:, :)
      !> array with eigenvalues
      complex(dp), intent(out)  :: omega(:)
      !> array with right eigenvectors
      complex(dp), intent(out)  :: vr(:, :)
    end subroutine qz_direct

    module subroutine arnoldi(matrix_A, matrix_B, omega, vr)
      !> matrix A
      complex(dp), intent(in)   :: matrix_A(:, :)
      !> matrix B
      real(dp), intent(in)      :: matrix_B(:, :)
      !> array with eigenvalues
      complex(dp), intent(out)  :: omega(:)
      !> array with right eigenvectors
      complex(dp), intent(out)  :: vr(:, :)
    end subroutine arnoldi
  end interface

  public  :: residual_norm
  public  :: solve_evp
  public  :: solvers_clean

contains


  !> Main subroutine to solve the eigenvalue problem. Depending on the solvelist
  !! passed in the parfile, different solvers are called.
  !! @warning Throws an error if an unknown solver is passed. @endwarning
  subroutine solve_evp(matrix_A, matrix_B, omega, vr)
    use mod_global_variables, only: solver

    !> A-matrix
    complex(dp), intent(in)   :: matrix_A(:, :)
    !> B-matrix
    real(dp), intent(in)      :: matrix_B(:, :)
    !> eigenvalues
    complex(dp), intent(out)  :: omega(:)
    !> right eigenvectors
    complex(dp), intent(out)  :: vr(:, :)

    select case(solver)
    case("QR-invert")
      call qr_invert(matrix_A, matrix_B, omega, vr)
    case("QZ-direct")
      call qz_direct(matrix_A, matrix_B, omega, vr)
    case("arnoldi")
      call arnoldi(matrix_A, matrix_B, omega, vr)
    case default
      call log_message("unknown solver passed: " // solver, level="error")
      return
    end select
  end subroutine solve_evp


  !> Cleanup routine.
  subroutine solvers_clean()
    if (allocated(residual_norm)) then
      deallocate(residual_norm)
    end if
  end subroutine solvers_clean

end module mod_solvers