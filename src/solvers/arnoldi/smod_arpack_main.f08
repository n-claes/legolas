! =============================================================================
!> Main module for the Arnoldi-type solvers. Contains interfaces to the general Arnoldi
!! procedures (general, shift-invert, etc.).
submodule (mod_solvers) smod_arpack_main
  use mod_arpack_type, only: arpack_t, new_arpack_config
  use mod_linear_systems, only: solve_linear_system_complex_banded, &
    solve_linear_system_complex_banded_LU, get_LU_factorisation_banded
  use mod_settings, only: settings_t
  implicit none

  interface
    !> Solves the eigenvalue problem using the Arnoldi general method.
    module subroutine solve_arpack_general( &
      arpack_cfg, matrix_A, matrix_B, settings, omega, vr &
    )
      !> arpack configuration
      type(arpack_t), intent(in) :: arpack_cfg
      !> matrix A
      type(matrix_t), intent(in) :: matrix_A
      !> matrix B
      type(matrix_t), intent(in) :: matrix_B
      !> settings object
      type(settings_t), intent(in) :: settings
      !> array with eigenvalues
      complex(dp), intent(out)  :: omega(:)
      !> array with right eigenvectors
      complex(dp), intent(out)  :: vr(:, :)
    end subroutine solve_arpack_general

    !> Solves the eigenvalue problem using the Arnoldi shift-invert method.
    module subroutine solve_arpack_shift_invert( &
      arpack_cfg, matrix_A, matrix_B, settings, omega, vr &
    )
      !> arpack configuration
      type(arpack_t), intent(in) :: arpack_cfg
      !> matrix A
      type(matrix_t), intent(in) :: matrix_A
      !> matrix B
      type(matrix_t), intent(in) :: matrix_B
      !> settings object
      type(settings_t), intent(in) :: settings
      !> array with eigenvalues
      complex(dp), intent(out) :: omega(:)
      !> array with right eigenvectors
      complex(dp), intent(out) :: vr(:, :)
    end subroutine solve_arpack_shift_invert
  end interface

contains


  module procedure arnoldi
#if _ARPACK_FOUND
    !> type containing parameters for arpack configuration
    type(arpack_t) :: arpack_cfg

    select case(settings%solvers%get_arpack_mode())
    case("general")
      call log_message("Arnoldi iteration, general mode", level="debug")
      arpack_cfg = new_arpack_config( &
        evpdim=matrix_A%matrix_dim, &
        mode=1, &
        bmat="I", &
        solver_settings=settings%solvers &
      )
      call solve_arpack_general(arpack_cfg, matrix_A, matrix_B, settings, omega, vr)
    case("shift-invert")
      call log_message("Arnoldi iteration, shift-invert mode", level="debug")
      arpack_cfg = new_arpack_config( &
        evpdim=matrix_A%matrix_dim, &
        mode=2, &
        bmat="I", &
        solver_settings=settings%solvers &
      )
      call solve_arpack_shift_invert( &
        arpack_cfg, matrix_A, matrix_B, settings, omega, vr &
      )
    case default
      call log_message( &
        "unknown mode for ARPACK: " // settings%solvers%get_arpack_mode(), &
        level="error" &
      )
      return
    end select

    call arpack_cfg%destroy()

#else
  call log_message( &
    "ARPACK was not found and/or CMake failed to link", level="warning" &
  )
  call log_message( &
    "unable to use 'arnoldi', try another solver!", level="error" &
  )
#endif
  end procedure arnoldi

end submodule smod_arpack_main
