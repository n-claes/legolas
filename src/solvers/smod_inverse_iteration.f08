! =============================================================================
!> Submodule containing the implementation of the inverse iteration algorithm.
!! TODO more docs
submodule (mod_solvers) smod_inverse_iteration
  use mod_global_variables, only: dim_subblock, maxiter, tolerance, sigma
  use mod_check_values, only: is_equal
  use mod_banded_matrix, only: banded_matrix_t, new_banded_matrix
  use mod_banded_matrix_hermitian, only: hermitian_banded_matrix_t
  use mod_transform_matrix, only: matrix_to_banded, matrix_to_hermitian_banded
  implicit none

contains

  !> Solves for one eigenvalue using inverse iteration.
  !! @warning Throws an error if <tt>matrix_A</tt> or <tt>matrix_B</tt>
  !!          is not a square matrix. @endwarning
  module procedure inverse_iteration
    !> sparse representation of B
    type(hermitian_banded_matrix_t) :: B_band
    !> sparse representation of A
    type(banded_matrix_t)           :: A_band
    !> LU decomposition of A - sigma * B
    type(banded_matrix_t)           :: LU
    !> pivoting array of the LU decomposition
    integer, allocatable          :: LU_ipiv(:)
    !> Number of sub/superdiagonals of A and B
    integer                       :: kd
    !> order of matrices A and B
    integer                       :: N
    !> right-eigenvector
    complex(dp), allocatable      :: x(:)
    !> iteration vectors
    complex(dp), allocatable      :: r(:)
    complex(dp), allocatable      :: s(:)
    !> current approximation
    complex(dp)                   :: ev
    !> info parameter
    integer                       :: info
    !> iteration count
    integer                       :: i, j
    !> flag to see if converged
    logical                       :: converged
    !> local variable for the tolerance
    real(dp)                      :: tol
    !> types of used BLAS functions
    real(dp)                      :: dznrm2
    complex(dp)                   :: zdotc
    integer                       :: idamax

    ! check input sanity
    if (.not. (matrix_A%matrix_dim == matrix_B%matrix_dim)) then
      call log_message("A or B not square, or not compatible", level="error")
    end if

    ! if maxiter is not set in the parfile it's still 0, default to 100
    if (maxiter == 0) then
      maxiter = 100
    else if (maxiter < 0) then
      call log_message( &
        "maxiter has to be positive, but is equal to " // str(maxiter), level="error" &
      )
      return
    end if

    if (is_equal(sigma, (0.0d0, 0.0d0))) then
      call log_message( &
        "inverse-iteration: sigma can not be equal to zero", &
        level="error" &
      )
      return
    end if ! LCOV_EXCL_STOP

    ! set array dimensions
    N = matrix_A%matrix_dim
    kd = 2*dim_subblock+1 ! at most 2 subblocks away from diag

    ! allocate iteration vectors
    allocate(x(N))
    allocate(r(N))
    allocate(s(N))

    ! get sparse version B
    call matrix_to_hermitian_banded(matrix=matrix_B, diags=kd, uplo="U", banded=B_band)
    ! get sparse version of A
    call matrix_to_banded(matrix=matrix_A, subdiags=kd, superdiags=kd, banded=A_band)

    ! get LU decompositon of A - sigma * B
    ! NOTE: the resulting LU decomposition has double the lower diagonals
    ! Compute A - sigma * B
    LU = new_banded_matrix(rows=N, cols=N, subdiags=2*kd, superdiags=kd)
    do j = 1, N
      do i = max(1, j - kd), min(N, j + kd)
        LU%AB(2*kd + 1 + i - j, j) = ( &
          A_band%get_element(i, j) - sigma * B_band%get_element(i, j) &
        )
      end do
    end do
    ! Compute LU decomposition of A - sigma * B
    allocate(LU_ipiv(N))
    call zgbtrf(N, N, kd, kd, LU%AB, LU%kl+LU%ku+1, LU_ipiv, info)
    if (info /= 0) then ! LCOV_EXCL_START
      call log_message( &
        "[A - sigma * B](" // str(info) // "," // str(info) // ") is zero", &
        level="warning" &
      )
    end if ! LCOV_EXCL_STOP

    ! start iteration
    tol = tolerance
    i = 0
    ev = sigma
    ! use solution of U x = 1 as initial guess
    x = 1.0_dp
    call ztbsv("U", "N", "N", N, kd+kd, LU%AB, 2*kd+kd+1, x, 1)

    converged = .false.
    do while (i <= maxiter .and. .not.converged)
      ! end of 'last' iteration
      ! r = B x
      call zhbmv( &
        B_band%uplo, &
        B_band%n, &
        B_band%kd, &
        (1.0_dp, 0.0_dp), &
        B_band%AB, &
        B_band%kd + 1, &
        x, 1, &
        (0.0_dp, 0.0_dp), &
        r, &
        1 &
      )
      ! s = A x
      call zgbmv( &
        "N", &
        A_band%m, &
        A_band%n, &
        A_band%kl, &
        A_band%ku, &
        (1.0_dp, 0.0_dp), &
        A_band%AB, &
        A_band%kl + A_band%ku + 1, &
        x, &
        1, &
        (0.0_dp, 0.0_dp), &
        s, &
        1 &
      )

      ! ev = x^H s / x^H r
      ev = zdotc(N, x, 1, s, 1) / zdotc(N, x, 1, r, 1)

      ! s = s - ev*r
      call zaxpy(N, -ev, r, 1, s, 1)

      ! reached low enough residual?
      if (dznrm2(N, s, 1) < abs(ev)*tol) then
        converged = .true.
        exit
      end if

      ! start of 'next' iteration
      i = i + 1

      ! solve (A - sigma B) x = r
      x = r
      call zgbtrs("N", LU%n, kd, kd, 1, LU%AB, 2*kd+kd+1, LU_ipiv, x, LU%n, info)

      ! normalize x
      call zdscal(N, 1.0_dp / dznrm2(N, x, 1), x, 1)
    end do

    ! if we did not converge, raise a warning
    if (.not.converged) then
      if (i == maxiter+1) then
        call log_message( &
          "Inverse iteration failed to converge! (maxiter reached)", level="warning" &
        )
        call log_message( &
          "number of iterations: " // str(maxiter), &
          level="warning", &
          use_prefix=.false. &
        )
      else
        call log_message( &
          "Inverse iteration failed to converge! (divergence)", level="warning" &
        )
        call log_message( &
          "number of iterations: " // str(i), &
          level="warning", &
          use_prefix=.false. &
        )
      end if
    end if

    ! write found eigenvalue
    omega = ev

    ! write eigenvector if requested
    if (should_compute_eigenvectors()) then
      ! make largest coefficient real
      i = idamax(N, abs(x), 1)
      call zscal(N, conjg(x(i))/abs(x(i)), x, 1)
      ! write eigenvector
      vr(:, 1) = x
    end if

    ! tear down
    deallocate(x)
    deallocate(r)
    deallocate(s)
    call B_band%destroy()
    call A_band%destroy()
    call LU%destroy()
    deallocate(LU_ipiv)
  end procedure inverse_iteration

end submodule smod_inverse_iteration
