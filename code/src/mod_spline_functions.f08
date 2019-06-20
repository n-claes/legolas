!
! MODULE: mod_spline_functions
!
!> @author
!> Niels Claes
!> niels.claes@kuleuven.be
!
!> @note
!> We don't use rj_low (= r_{j-1}), rj_center (= r_j) or rj_high (= r_{j+1}),
!! because in every grid block [i, i+1] there is a contribution from both
!! the current functions in this interval as from the functions in the grid
!! block [i-1, i]. See page 193 of advanded MHD by Keppens and Poedts, this
!! holds both for cubic and quadratic functions.
!! Hence we use rj_lo for r_j and rj_hi for r_{j+1}, which in the interval
!! [i, i+1] (2 and 4) corresponds to the current quadratic/cubic functions,
!! while in the interval [i-1, i] (1 and 3) this corresponds to the previous
!! quadratic/cubic functions.
!
! DESCRIPTION:
!> Module to calculate the quadratic and cubic spline functions and their
!! derivatives.
!
module mod_spline_functions
  use mod_global_variables
  implicit none

  public

contains

  !> Calculates the quadratic basis functions.
  !! @param[in] r       The current value for r in grid_gauss
  !! @param[in] rj_lo   Left value of r in grid
  !! @param[in] rj_hi   Right value of r in grid (so r has 4 values between
  !!                    rj_hi and rj_lo)
  !! @param[out] h_quadratic  Array containing the 4 quadratic basis functions
  subroutine quadratic_factors(r, rj_lo, rj_hi, h_quadratic)
    real(dp), intent(in)  ::  r, rj_lo, rj_hi
    real(dp), intent(out) ::  h_quadratic(4)

    h_quadratic(1) = 4.0d0 * (r - rj_lo) * (rj_hi - r) / (rj_hi - rj_lo)**2
    h_quadratic(2) = 0.0d0
    h_quadratic(3) = (2.0d0*r - rj_hi - rj_lo) * (r - rj_lo) / (rj_hi - rj_lo)**2
    h_quadratic(4) = (2.0d0*r - rj_hi - rj_lo) * (r - rj_hi) / (rj_hi - rj_lo)**2

  end subroutine quadratic_factors

  !> Calculates the derivatives of the quadratic basis functions.
  !! @param[in] r       The current value for r in grid_gauss
  !! @param[in] rj_lo   Left value of r in grid
  !! @param[in] rj_hi   Right value of r in grid (so r has 4 values between
  !!                    rj_hi and rj_lo)
  !! @param[out] dh_quadratic_dr  Array containing the derivatives of the
  !!                              4 quadratic basis functions
  subroutine quadratic_factors_deriv(r, rj_lo, rj_hi, dh_quadratic_dr)
    real(dp), intent(in)  ::  r, rj_lo, rj_hi
    real(dp), intent(out) ::  dh_quadratic_dr(4)

    dh_quadratic_dr(1) = 4.0d0 * (-2.0d0*r + rj_hi + rj_lo) / (rj_hi - rj_lo)**2
    dh_quadratic_dr(2) = 0.0d0
    dh_quadratic_dr(3) = (4.0d0*r - rj_hi - 3.0d0*rj_lo) / (rj_hi - rj_lo)**2
    dh_quadratic_dr(4) = (4.0d0*r - rj_lo - 3.0d0*rj_hi) / (rj_hi - rj_lo)**2

  end subroutine quadratic_factors_deriv

  !> Calculates the cubic basis functions.
  !! @param[in] r       The current value for r in grid_gauss
  !! @param[in] rj_lo   Left value of r in grid
  !! @param[in] rj_hi   Right value of r in grid (so r has 4 values between
  !!                    rj_hi and rj_lo)
  !! @param[out] h_cubic  Array containing the 4 cubic basis functions
  subroutine cubic_factors(r, rj_lo, rj_hi, h_cubic)
    real(dp), intent(in)  :: r, rj_lo, rj_hi
    real(dp), intent(out) :: h_cubic(4)

    h_cubic(1) =  3.0d0 * ( (r - rj_lo) / (rj_hi - rj_lo) )**2 &
                 -2.0d0 * ( (r - rj_lo) / (rj_hi - rj_lo) )**3
    h_cubic(2) =  3.0d0 * ( (rj_hi - r) / (rj_hi - rj_lo) )**2 &
                 -2.0d0 * ( (rj_hi - r) / (rj_hi - rj_lo) )**3
    h_cubic(3) = (r - rj_hi) * ( (r - rj_lo) / (rj_hi - rj_lo) )**2
    h_cubic(4) = (r - rj_lo) * ( (rj_hi - r) / (rj_hi - rj_lo) )**2

  end subroutine cubic_factors

  !> Calculates the derivatives of the cubic basis functions.
  !! @param[in] r       The current value for r in grid_gauss
  !! @param[in] rj_lo   Left value of r in grid
  !! @param[in] rj_hi   Right value of r in grid (so r has 4 values between
  !!                    rj_hi and rj_lo)
  !! @param[out] dh_cubic_dr  Array containing the derivatives of the
  !!                          4 cubic basis functions
  subroutine cubic_factors_deriv(r, rj_lo, rj_hi, dh_cubic_dr)
    real(dp), intent(in)  :: r, rj_lo, rj_hi
    real(dp), intent(out) :: dh_cubic_dr(4)

    dh_cubic_dr(1) =  6.0d0 * (r - rj_lo) / (rj_hi - rj_lo)**2 &
                     -6.0d0 * (r - rj_lo)**2 / (rj_hi - rj_lo)**3
    dh_cubic_dr(2) = -6.0d0 * (rj_hi - r) / (rj_hi - rj_lo)**2 &
                     +6.0d0 * (rj_hi - r)**2 / (rj_hi - rj_lo)**3
    dh_cubic_dr(3) = ( 2.0d0*(r - rj_hi) * (r - rj_lo) + (r - rj_lo)**2 ) &
                     / (rj_hi - rj_lo)**2
    dh_cubic_dr(4) = ( 2.0d0*(r - rj_lo) * (r - rj_hi) + (r - rj_hi)**2 ) &
                     / (rj_hi - rj_lo)**2

  end subroutine cubic_factors_deriv





end module mod_spline_functions
