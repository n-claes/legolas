module mod_setup_matrix_b
  implicit none

  ! Sets up the B-matrix for the eigenvalue problem wBX = AX
  real                      :: h_cubic(4), h_quadratic(4)

contains

  subroutine construct_B(matrix_gridpts, matrix_B)
    integer, intent(in)     ::  matrix_gridpts
    real, intent(inout)     ::  matrix_B(matrix_gridpts, matrix_gridpts)

    integer                 :: i, j

    ! Initialize matrix to zero
    do i=1, matrix_gridpts
      do j=1, matrix_gridpts
        matrix_B(i, j) = 0.0d0
      end do
    end do

    print*,matrix_B


  end subroutine construct_B



end module mod_setup_matrix_b
