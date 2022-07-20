! =============================================================================
!> Module containing the implementation for the ARPACK general-type solver, that is,
!! given the general eigenvalue problem $$ AX = \omega BX, $$ find \(k\) eigenvalues
!! that satisfy a given criterion.
submodule (mod_solvers:smod_arpack_main) smod_arpack_general
  use mod_banded_matrix, only: banded_matrix_t
  use mod_transform_matrix, only: matrix_to_banded
  implicit none

contains

  !> Implementation of the ARPACK general solver.
  module procedure solve_arpack_general
    !> contains the basis vectors
    complex(dp) :: basis_vectors(arpack_cfg%get_evpdim(), arpack_cfg%get_ncv())
    !> work array of length 3N
    complex(dp) :: workd(3 * arpack_cfg%get_evpdim())
    !> work array
    complex(dp) :: workl(arpack_cfg%get_lworkl())
    !> work array of length ncv
    real(dp) :: rwork(arpack_cfg%get_ncv())
    !> integer array with pointers to mark work array locations
    integer :: ipntr(14)
    !> logical array of dimension ncv, sets which Ritz vectors to compute
    logical :: select_vectors(arpack_cfg%get_ncv())
    !> work array of length 2*ncv
    complex(dp) :: workev(2 * arpack_cfg%get_ncv())

    integer :: diags
    logical :: converged
    type(banded_matrix_t) :: bmat_banded
    integer :: xstart, xend, ystart, yend

    ! we fill at most the full quadblock, so -1 diagonals
    ! IDEA: maybe do a pass over the B-matrix first to figure out how many diagonals
    ! we need? This is problem-dependent so in some cases there might be quite some
    ! room here for optimisations
    diags = dim_quadblock - 1
    call log_message("converting B-matrix into banded structure", level="debug")
    call matrix_to_banded( &
      matrix=matrix_B, subdiags=diags, superdiags=diags, banded=bmat_banded &
    )

    call log_message("doing Arnoldi iteration", level="debug")
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

      ! in the following y is given by workd(ipntr(2)), x by workd(ipntr(1)).
      ! we need to explicitly set start and end to avoid incompatible ranks
      xstart = ipntr(1)
      xend = xstart + arpack_cfg%get_evpdim() - 1
      ystart = ipntr(2)
      yend = ystart + arpack_cfg%get_evpdim() - 1

      ! note that the linked-list datastructures have the matrix-vector product
      ! implemented using the operator (*).
      select case(arpack_cfg%ido)
      case(-1, 1)
        ! get y <--- OP*x, we do not calculate OP*x explicitly.
        ! We need R = OP*x = inv[B]*A*x, so do the following:
        ! 1. calculate u = A*x
        ! 2. solve linear system B * R = u for R
        workd(ystart:yend) = solve_linear_system_complex_banded( &
          bandmatrix=bmat_banded, vector=matrix_A * workd(xstart:xend) &
        )
      case default
        ! when convergence is achieved or maxiter is reached
        exit
      end select
    end do

    call arpack_cfg%parse_znaupd_info(converged)

    ! if we have a normal exit, extract the eigenvalues through zneupd
    call zneupd( &
      .true., &  ! always calculate eigenvectors, negligible additional cost in ARPACK
      "A", &  ! calculate Ritz vectors
      select_vectors, &
      omega(1:arpack_cfg%get_nev()), &
      vr(:, 1:arpack_cfg%get_nev()), &
      size(vr, dim=1), &
      (0.0d0, 0.0d0), &  ! sigma value, not needed here
      workev, &
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

    call arpack_cfg%parse_zneupd_info()
    call arpack_cfg%parse_finished_stats()
  end procedure solve_arpack_general

end submodule smod_arpack_general
