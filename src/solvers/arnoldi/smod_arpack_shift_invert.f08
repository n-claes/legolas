submodule (mod_solvers:smod_arpack_main) smod_arpack_shift_invert
  use mod_banded_matrix, only: banded_matrix_t
  use mod_transform_matrix, only: matrix_to_banded, banded_to_array
  implicit none

contains

  module procedure solve_arpack_shift_invert
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
    type(banded_matrix_t) :: amat_min_sigmab
    integer :: xstart, xend, ystart, yend

    diags = dim_quadblock - 1
    call log_message("converting A - sigma*B into banded structure", level="debug")
    call matrix_to_banded( &
      matrix=matrix_A - matrix_B * sigma, &
      subdiags=diags, &
      superdiags=diags, &
      banded=amat_min_sigmab &
    )

    call log_message("doing Arnoldi shift-invert", level="debug")
    converged = .false.

    do while(.not. converged)
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

      ! x is given by workd(ipntr(1))
      xstart = ipntr(1)
      xend = xstart + arpack_cfg%get_evpdim() - 1
      ! y is given by workd(ipntr(2))
      ystart = ipntr(2)
      yend = ystart + arpack_cfg%get_evpdim() - 1

      select case(arpack_cfg%ido)
      case(-1, 1)
        ! ido = -1 on first call, forces starting vector in OP range
        ! get y <--- OP * x
        ! we need R = OP*x = inv[A - sigma*B]*B*x
        ! 1. calculate u = B*x
        ! 2. solve linear system [A - sigma*B] * R = u for R
        workd(ystart:yend) = solve_linear_system_complex_banded( &
          bandmatrix=amat_min_sigmab, vector=matrix_B * workd(xstart:xend) &
        )
      case default
        ! when convergence is achieved or maxiter is reached
        exit
      end select
    end do

    ! check info parameter from znaupd, this errors if necessary
    call arpack_cfg%parse_znaupd_info(converged)
    ! if we have a normal exit, extract the eigenvalues through zneupd
    call zneupd( &
      .true., &  ! always calculate eigenvectors, negligible additional cost in ARPACK
      "A", &  ! calculate Ritz vectors
      select_vectors, &
      omega(1:arpack_cfg%get_nev()), &
      vr(:, 1:arpack_cfg%get_nev()), &
      size(vr, dim=1), &
      sigma, &
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
    !> @note In applying shift-invert we made the transformation $C = inv[B]*A$ and
    !! solved the standard eigenvalue problem $Cx = \nu x$ instead since B isn't
    !! always Hermitian (e.g. if we include Hall).
    !! According to the ARPACK documentation, section 3.2.2, this
    !! implies that we must manually transform the eigenvalues $\nu_j$ from $C$ to the
    !! eigenvalues $\omega_j$ from the original system. This uses the relation
    !! $$ \omega_j = \sigma + \frac{1}{\nu_j} $$
    !! @endnote
    omega = sigma + (1.0d0 / omega)

    call arpack_cfg%parse_zneupd_info()
    call arpack_cfg%parse_finished_stats()
  end procedure solve_arpack_shift_invert

end submodule smod_arpack_shift_invert
