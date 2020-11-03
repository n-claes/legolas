module mod_solvers
  use mod_global_variables, only: dp, matrix_gridpts, write_eigenfunctions
  use mod_logging, only: log_message, char_log, dp_fmt, int_fmt
  use mod_check_values, only: check_small_values
  implicit none

  private

  !> interface to the different solution methods implemented in submodules
  interface
    module subroutine qr_invert(matrix_A, matrix_B, omega, vl, vr)
      complex(dp), intent(in)   :: matrix_A(:, :)
      real(dp), intent(in)      :: matrix_B(:, :)
      complex(dp), intent(out)  :: omega(:)
      complex(dp), intent(out)  :: vl(:, :)
      complex(dp), intent(out)  :: vr(:, :)
    end subroutine

    module subroutine qz_direct(matrix_A, matrix_B, omega, vl, vr)
      complex(dp), intent(in)   :: matrix_A(:, :)
      real(dp), intent(in)      :: matrix_B(:, :)
      complex(dp), intent(out)  :: omega(:)
      complex(dp), intent(out)  :: vl(:, :)
      complex(dp), intent(out)  :: vr(:, :)
    end subroutine qz_direct
  end interface

  public  :: solve_evp

contains


  subroutine solve_evp(matrix_A, matrix_B, omega, vl, vr)
    use mod_global_variables, only: solver

    complex(dp), intent(in)   :: matrix_A(:, :)
    real(dp), intent(in)      :: matrix_B(:, :)
    complex(dp), intent(out)  :: omega(:)
    complex(dp), intent(out)  :: vl(:, :)
    complex(dp), intent(out)  :: vr(:, :)

    select case(solver)
      case("QR-invert")
        call qr_invert(matrix_A, matrix_B, omega, vl, vr)
      case("QZ-direct")
        call qz_direct(matrix_A, matrix_B, omega, vl, vr)
      case default
        call log_message("unknown solver passed: " // solver, level="error")
    end select

  end subroutine solve_evp

end module mod_solvers