module mod_banded_operations
  use mod_global_variables, only: dp
  use mod_banded_matrix, only: banded_matrix_t
  implicit none

  private

  interface multiply
    module procedure banded_matrix_x_vector
  end interface multiply

  public :: multiply

contains

  !> Calculates the matrix-vector product of a general complex banded matrix and a
  !! complex vector. Uses the level 2 BLAS routine <tt>zgbmv</tt>.
  function banded_matrix_x_vector(bandmatrix, vector) result(rvector)
    !> the banded matrix
    type(banded_matrix_t), intent(in) :: bandmatrix
    !> the vector
    complex(dp), intent(in) :: vector(bandmatrix%n)
    !> the resulting matrix-vector product
    complex(dp) :: rvector(bandmatrix%m)

    call zgbmv( &
      "N", &
      bandmatrix%m, &
      bandmatrix%n, &
      bandmatrix%kl, &
      bandmatrix%ku, &
      (1.0_dp, 0.0_dp), &
      bandmatrix%AB, &
      size(bandmatrix%AB, dim=1), &
      vector, &
      1, &
      (0.0_dp, 0.0_dp), &
      rvector, &
      1 &
    )
  end function banded_matrix_x_vector

end module mod_banded_operations
