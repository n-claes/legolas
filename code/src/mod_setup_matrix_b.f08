module mod_setup_matrix_b
  implicit none

  ! Sets up the B-matrix for the eigenvalue problem wBX = AX
  real                      :: h_cubic(4), h_quadratic(4)

contains

  subroutine construct_B(grid, matrix_B)
    use mod_global_variables

    real, intent(in)        :: grid(integral_gridpts)
    real, intent(inout)     :: matrix_B(matrix_gridpts, matrix_gridpts)
    integer                 :: i, j
    double precision        :: r_lo, r, r_hi, eps, d_eps_dr

    ! Initialize matrix to zero
    do i=1, matrix_gridpts
      do j=1, matrix_gridpts
        matrix_B(i, j) = 0.0d0
      end do
    end do

    do i = 2, integral_gridpts-1
      r_lo = grid(i - 1)
      r    = grid(i)
      r_hi = grid(i + 1)

      if (geometry .eq. "cartesian") then
        eps      = 1.0d0
        d_eps_dr = 0.0d0
      else if (geometry .eq. "cylindrical") then
        !TODO: IN LEDA EPS EQUALS 'ZSR' (WHICH SHOULD BE R)
        !      FIGURE OUT WHAT THIS SHOULD BE BEFORE GOING FURTHER
        eps      = 0
        d_eps_dr = 1.0d0
      endif
    end do

  end subroutine construct_B



end module mod_setup_matrix_b
