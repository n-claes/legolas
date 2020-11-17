module mod_solvers
  use mod_global_variables, only: dp, matrix_gridpts, write_eigenfunctions
  use mod_logging, only: log_message, char_log, char_log2, dp_fmt, int_fmt
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

    module subroutine arnoldi(matrix_A, matrix_B, omega, vr)
      complex(dp), intent(in)   :: matrix_A(:, :)
      real(dp), intent(in)      :: matrix_B(:, :)
      complex(dp), intent(out)  :: omega(:)
      complex(dp), intent(out)  :: vr(:, :)
    end subroutine arnoldi
  end interface

  public  :: solve_evp
  public  :: do_arpack_sanity_checks

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
    case("arnoldi")
      call arnoldi(matrix_A, matrix_B, omega, vr)
    case default
      call log_message("unknown solver passed: " // solver, level="error")
      return
    end select

  end subroutine solve_evp


  subroutine do_arpack_sanity_checks(evpdim)
    use mod_global_variables, only: number_of_eigenvalues, which_eigenvalues, maxiter

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

end module mod_solvers