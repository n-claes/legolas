submodule (mod_solvers:smod_arpack_main) smod_arpack_general
  implicit none

contains

  module procedure solve_arpack_general
    complex(dp) :: basis_vectors(arpack_cfg%get_evpdim(), arpack_cfg%get_ncv())
    complex(dp) :: workd(3 * arpack_cfg%get_evpdim())
    complex(dp) :: workl(arpack_cfg%get_lworkl())
    real(dp) :: rwork(arpack_cfg%get_ncv())
    integer :: ipntr(14)
    logical :: converged

    return

    call log_message("starting Arnoldi iteration", level="debug")

    converged = .false.
    ! we keep iterating for as long as the eigenvalues are not converged.
    ! If convergence is achieved or the maximum number of iterations is reached,
    ! ARPACK sets ido=99 so the while-loop will always break.
    do while (.not. converged)
      call znaupd( &
        arpack_cfg%ido, &
        arpack_cfg%get_bmat(), &
        arpack_cfg%get_evpdim(), &
        arpack_cfg%get_which(), &
        arpack_cfg%get_nev(), &
        arpack_cfg%get_tolerance(), &
        arpack_cfg%residual, &
        arpack_cfg%get_ncv(), &
        basis_vectors, &
        size(basis_vectors, dim=1), &
        arpack_cfg%iparam, &
        ipntr, &
        workd, &
        workl, &
        arpack_cfg%get_lworkl(), &
        rwork, &
        arpack_cfg%info &
      )
    end do
  end procedure solve_arpack_general

end submodule smod_arpack_general