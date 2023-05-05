! =============================================================================
!> Module containing functions to solve linear systems $$ AX = B $$.
module mod_linear_systems
  use mod_global_variables, only: dp
  use mod_banded_matrix, only: banded_matrix_t
  use mod_logging, only: logger, str
  implicit none

  private

  public :: solve_linear_system_complex_banded
  public :: solve_linear_system_complex_banded_LU
  public :: get_LU_factorisation_banded

contains

  !> Calculates the solution \(X\) to a system of linear equations \(AX = B\) where
  !! \(A\) is a complex banded matrix and \(B\) is a complex vector.
  function solve_linear_system_complex_banded(bandmatrix, vector) result(xvector)
    !> the banded \(A\)-matrix on the left-hand side
    type(banded_matrix_t), intent(in) :: bandmatrix
    !> the \(B\)-vector on the right-hand side
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
    do j = 1, bandmatrix%n
      do i = max(1, j - ku), min(nb_eqs, j + kl)
        ABmat(kl + ku + 1 + i - j, j) = bandmatrix%get_element(row=i, col=j)
      end do
    end do
    ! on entry, xvector equals the B-vector
    xvector = vector
    ! on exit, xvector contains solution if info == 0
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
    if (info /= 0) call throw_info_nonzero_warning(info, "zgbsv")
    deallocate(ABmat)
  end function solve_linear_system_complex_banded


  !> Calculates the solution \(X\) to a system of linear equations \(AX = B\) where
  !! \(A\) is a complex banded matrix and \(B\) is a complex vector.
  !! Uses the LU factorisation of \(A\) and pivoting information from zgbtrf.
  function solve_linear_system_complex_banded_LU(bandmatrix, vector, LU, ipiv) &
    result(xvector)
    !> the banded \(A\)-matrix on the left-hand side
    type(banded_matrix_t), intent(in) :: bandmatrix
    !> the \(B\)-vector on the right-hand side
    complex(dp), intent(in) :: vector(:)
    !> the LU factorisation of \(A\), from zgbtrf
    complex(dp), intent(in) :: LU(:, :)
    !> the pivoting information from zgbtrf
    integer, intent(in) :: ipiv(:)
    !> the solution vector \(X\)
    complex(dp) :: xvector(size(vector))
    integer :: info

    ! on entry, xvector equals the B-vector
    xvector = vector
    call zgbtrs( &
      "N", &
      bandmatrix%n, &
      bandmatrix%kl, &
      bandmatrix%ku, &
      1, &
      LU, &
      size(LU, dim=1), &
      ipiv, &
      xvector, &
      size(vector), &
      info &
    )
    if (info /= 0) call throw_info_nonzero_warning(info, "zgbtrs")
  end function solve_linear_system_complex_banded_LU


  !> Calculates the LU factorisation of a complex banded matrix \(A\).
  !! Uses the LAPACK routine zgbtrf.
  subroutine get_LU_factorisation_banded(bandmatrix, LU, ipiv)
    !> the banded \(A\)-matrix on the left-hand side
    type(banded_matrix_t), intent(in) :: bandmatrix
    !> the LU factorisation of \(A\)
    complex(dp), allocatable, intent(out) :: LU(:, :)
    !> the pivoting information after factorisation
    integer, allocatable, intent(out) :: ipiv(:)
    integer :: nbrows, nbcols, kl, ku, irow, icol, info

    nbrows = bandmatrix%m
    nbcols = bandmatrix%n
    kl = bandmatrix%kl
    ku = bandmatrix%ku
    ! Fill LU using the bandmatrix, A(i, j) is mapped to LU(kl + ku + 1 + i - j, j)
    ! for max(1, j - ku) <= i <= min(nbrows, j + kl)
    allocate(LU(2 * kl + ku + 1, nbrows))
    do icol = 1, nbcols
      do irow = max(1, icol - ku), min(nbrows, icol + kl)
        LU(kl + ku + 1 + irow - icol, icol) = bandmatrix%get_element(row=irow, col=icol)
      end do
    end do
    allocate(ipiv(min(nbrows, nbcols)))

    call zgbtrf(nbrows, nbcols, kl, ku, LU, size(LU, dim=1), ipiv, info)
    if (info /= 0) call throw_info_nonzero_warning(info, "zgbtrf")
  end subroutine get_LU_factorisation_banded


  ! LCOV_EXCL_START
  subroutine throw_info_nonzero_warning(info, routine_name)
    integer, intent(in) :: info
    character(len=*), intent(in) :: routine_name

    call logger%warning( &
      "LAPACK routine " // routine_name // " failed! info = " // str(info) &
    )
  end subroutine throw_info_nonzero_warning
  ! LCOV_EXCL_STOP

end module mod_linear_systems
