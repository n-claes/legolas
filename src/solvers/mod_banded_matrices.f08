! =============================================================================
!> Module implementing a type and utilities for LAPACK band storage matrices.
module mod_banded_matrices
#include <assert.fpp>
  use mod_global_variables, only: dp
  implicit none

  !> A normal complex banded matrix
  type banded_matrix
    !> Number of rows
    integer     :: m
    !> Number of columns
    integer     :: n
    !> Number of subdiagonals
    integer     :: kl
    !> Number of superdiagonals
    integer     :: ku
    !> The matrix in banded storage
    complex(dp), allocatable :: BS(:, :)
  end type banded_matrix

  !> A Hermitain complex banded matrix
  type hermitian_banded_matrix
    !> Number of rows/columns
    integer     :: n
    !> If 'U' the upper triangular part is stored,
    !! if 'L' the lower triangular part is stored
    character   :: uplo
    !> Number of sub/superdiagonals
    integer     :: kd
    !> The matrix in banded storage
    complex(dp), allocatable :: BS(:, :)
  end type hermitian_banded_matrix

  private

  public :: banded_matrix
  public :: allocate_banded_matrix
  public :: deallocate_banded_matrix
  public :: real_dense_to_banded
  public :: complex_dense_to_banded
  public :: banded_to_complex_dense

  public :: hermitian_banded_matrix
  public :: allocate_hermitian_banded_matrix
  public :: deallocate_hermitian_banded_matrix
  public :: real_dense_to_hermitian_banded
  public :: complex_dense_to_hermitian_banded

contains

  !> Allocate a general banded matrix, frees existing
  !! memory if needed.
  subroutine allocate_banded_matrix(m, n, kl, ku, AB)
    !> Number of rows
    integer, intent(in) :: m
    !> Number of columns
    integer, intent(in) :: n
    !> Number subdiagonals
    integer, intent(in) :: kl
    !> Number superdiagonals
    integer, intent(in) :: ku
    !> Allocated banded matrix
    type(banded_matrix), intent(out) :: AB

    AB%m  = m;  AB%n  = n
    AB%kl = kl; AB%ku = ku
    if (allocated(AB%BS)) then
      deallocate(AB%BS)
    end if
    allocate(AB%BS(kl + ku + 1, n))
  end subroutine allocate_banded_matrix

  !> Deallocate a general banded matrix, if allocated
  subroutine deallocate_banded_matrix(AB)
    !> The banded matrix to deallocate
    type(banded_matrix), intent(inout) :: AB

    if (allocated(AB%BS)) then
      deallocate(AB%BS)
    end if
  end subroutine deallocate_banded_matrix

  !> Extract complex general banded matrix from a
  !! real dense matrix, allocating output if needed.
  subroutine real_dense_to_banded(A, kl, ku, AB)
    !> The source matrix
    real(dp), intent(in) :: A(:, :)
    !> Number of subdiagonals to extract
    integer, intent(in)  :: kl
    !> Number of superdiagonals to extract
    integer, intent(in)  :: ku
    !> The resulting complex banded matrix
    type(banded_matrix), intent(inout) :: AB

    integer :: m, n, i, j

    m = size(A, 1)
    n = size(A, 2)
    call allocate_banded_matrix(m, n, kl, ku, AB)
    do j = 1, n
      do i = max(1, j - ku), min(m, j + kl)
        AB%BS(ku + 1 + i - j, j) = cmplx(A(i, j), kind=dp)
      end do
    end do
  end subroutine real_dense_to_banded

  !> Extract complex general banded matrix from a
  !! complex dense matrix, allocating output if needed.
  subroutine complex_dense_to_banded(A, kl, ku, AB)
    !> The source matrix
    complex(dp), intent(in) :: A(:, :)
    !> Number of subdiagonals to extract
    integer, intent(in)     :: kl
    !> Number of superdiagonals to extract
    integer, intent(in)     :: ku
    !> The resulting complex banded matrix
    type(banded_matrix), intent(inout) :: AB

    integer :: m, n, i, j

    m = size(A, 1)
    n = size(A, 2)
    call allocate_banded_matrix(m, n, kl, ku, AB)
    do j = 1, n
      do i = max(1, j - ku), min(m, j + kl)
        AB%BS(ku + 1 + i - j, j) = A(i, j)
      end do
    end do
  end subroutine complex_dense_to_banded

  !> Expand complex general banded matrix to a
  !! dense complex matrix.
  subroutine banded_to_complex_dense(AB, A)
    !> The source banded matrix
    type(banded_matrix), intent(in) :: AB
    !> The resulting dense matrix
    complex(dp), intent(out) :: A(AB%m, AB%n)

    integer :: m, n, kl, ku, i, j

    m  = AB%m;  n  = AB%n
    kl = AB%kl; ku = AB%ku
    do j = 1, n
      do i = 1, m
        if (max(1, j - ku) <= i .and. &
            i <= min(m, j + kl)) then
          A(i, j) = AB%BS(ku + 1 + i - j, j)
        else
          A(i, j) = (0.0_dp, 0.0_dp)
        end if
      end do
    end do
  end subroutine banded_to_complex_dense


  !> Allocate a Hermitian banded matrix, frees existing
  !! memory if needed.
  subroutine allocate_hermitian_banded_matrix(n, uplo, kd, AB)
    !> Number of rows/columns
    integer, intent(in)   :: n
    !> If 'U' the upper triangular part is stored,
    !! if 'L' the lower triangular part is stored
    character, intent(in) :: uplo
    !> Number of sub/superdiagonals
    integer, intent(in)   :: kd
    !> Allocated banded Hermitian matrix
    type(hermitian_banded_matrix), intent(out) :: AB

    AB%n    = n
    AB%uplo = uplo
    AB%kd   = kd
    if (allocated(AB%BS)) then
      deallocate(AB%BS)
    end if
    allocate(AB%BS(kd + 1, n))
  end subroutine allocate_hermitian_banded_matrix

  !> Deallocate a Hermitian banded matrix, if allocated
  subroutine deallocate_hermitian_banded_matrix(AB)
    !> The Hermitian banded matrix to deallocate
    type(hermitian_banded_matrix), intent(inout) :: AB

    if (allocated(AB%BS)) then
      deallocate(AB%BS)
    end if
  end subroutine deallocate_hermitian_banded_matrix

  !> Extract Hermitian banded matrix from a real dense matrix,
  !! allocating output if needed.
  subroutine real_dense_to_hermitian_banded(A, uplo, kd, AB)
    !> The source matrix
    real(dp), intent(in)  :: A(:, :)
    !> If 'U' the upper triangular part is extracted,
    !! if 'L' the lower triangular part is extracted
    character, intent(in) :: uplo
    !> Number of sub/superdiagonals to extract
    integer, intent(in)   :: kd
    !> The resulting banded Hermitian matrix
    type(hermitian_banded_matrix), intent(inout) :: AB

    integer :: n, i, j

    assert(size(A, 1) == size(A, 2))
    assert(uplo == 'U' .or. uplo == 'L')

    n = size(A, 1)
    call allocate_hermitian_banded_matrix(n, uplo, kd, AB)

    if (uplo == 'U') then
      do j = 1, n
        do i = max(1, j - kd), j
          AB%BS(kd + 1 + i - j, j) = cmplx(A(i, j), kind=dp)
        end do
      end do
    else ! uplo == 'L'
      do j = 1, n
        do i = j, min(n, j + kd)
          AB%BS(1 + i - j, j) = A(i, j)
        end do
      end do
    end if
  end subroutine real_dense_to_hermitian_banded

  !> Extract Hermitian banded matrix from a complex dense matrix,
  !! allocating output if needed.
  subroutine complex_dense_to_hermitian_banded(A, uplo, kd, AB)
    !> The source matrix
    complex(dp), intent(in) :: A(:, :)
    !> If 'U' the upper triangular part is extracted,
    !! if 'L' the lower triangular part is extracted
    character, intent(in)   :: uplo
    !> Number of sub/superdiagonals to extract
    integer, intent(in)     :: kd
    !> The resulting banded Hermitian matrix
    type(hermitian_banded_matrix), intent(inout) :: AB

    integer :: n, i, j

    assert(size(A, 1) == size(A, 2))
    assert(uplo == 'U' .or. uplo == 'L')

    n = size(A, 1)
    call allocate_hermitian_banded_matrix(n, uplo, kd, AB)

    if (uplo == 'U') then
      do j = 1, n
        do i = max(1, j - kd), j
          AB%BS(kd + 1 + i - j, j) = A(i, j)
        end do
      end do
    else ! uplo == 'L'
      do j = 1, n
        do i = j, min(n, j + kd)
          AB%BS(1 + i - j, j) = A(i, j)
        end do
      end do
    end if
  end subroutine complex_dense_to_hermitian_banded

end module mod_banded_matrices
