! =============================================================================
!> Module containing functions to solve linear systems $$ AX = B $$.
module mod_linear_systems
  use mod_global_variables, only: dp
  use mod_banded_matrix, only: banded_matrix_t
  use mod_logging, only: log_message, str
  implicit none

  private

  public :: solve_linear_system_complex_banded

contains

  !> Calculates the solution \(X\) to a system of linear equations \(AX = B\) where
  !! \(A\) is a complex banded matrix and \(B\) is a complex vector.
  function solve_linear_system_complex_banded(bandmatrix, vector) result(xvector)
    !> the banded \(A\)-matrix on the left-hand side
    type(banded_matrix_t), intent(in) :: bandmatrix
    !> the \(B\}-vector on the right-hand side
    complex(dp), intent(in) :: vector(:)
    !> the solution vector \(X\)
    complex(dp) :: xvector(size(vector))
    complex(dp), allocatable :: ABmat(:, :)
    integer :: nb_eqs, kl, ku, i, j, info
    integer :: ipiv(bandmatrix%m)

    ! number of linear equations
    nb_eqs = bandmatrix%m
    ! number of subdiagonals
    kl = bandmatrix%kl
    ! number of superdiagonals
    ku = bandmatrix%ku
    ! ABmat needs an additional kl superdiagonals for internal row interchanges.
    ! A(i, j) is mapped to ABmat(kl + ku + 1 + i - j, j)
    allocate(ABmat(2 * kl + ku + 1, nb_eqs))
    do j = 1, bandmatrix%m
      do i = max(1, j - ku), min(nb_eqs, j + kl)
        ABmat(kl + ku + 1 + i - j, j) = bandmatrix%get_element(row=i, col=j)
      end do
    end do
    ! on entry, xvector equals the B-vector
    xvector = vector

    ! on exit, xvector contains solution (if info == 0)
    call zgbsv( &
      nb_eqs, &
      kl, &
      ku, &
      1, &
      ABmat, &
      size(ABmat, dim=1), &
      ipiv, &
      xvector, &
      size(vector), &
      info &
    )

    if (info /= 0) then  ! LCOV_EXCL_START
      call log_message( &
        "LAPACK routine zgbsv failed! info = " // str(info), level="warning" &
      )
    end if ! LCOV_EXCL_STOP
    deallocate(ABmat)
  end function solve_linear_system_complex_banded

end module mod_linear_systems
