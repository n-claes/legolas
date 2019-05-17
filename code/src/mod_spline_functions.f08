module mod_spline_functions
  implicit none

  ! spline functions and derivatives for quadratic and cubic elements
  ! rj_lo equals r_{j - 1}, rj_hi equals r_{j + 1}


contains

  subroutine quadratic_factors(r, rj_lo, rj, rj_hi, h_quadratic)
    integer, intent(in)           ::  r, rj_lo, rj, rj_hi
    real, intent(out)             ::  h_quadratic(4)

    h_quadratic(1) = 4.0d0 * (r - rj_lo) * (rj - r) / (rj - rj_lo)**2
    h_quadratic(2) = 0.0d0
    h_quadratic(3) = (2.0d0*r - rj - rj_lo) * (r - rj_lo) / (rj - rj_lo)**2
    h_quadratic(4) = (2.0d0*r - rj_hi - rj) * (r - rj_hi) / (rj_hi - rj)**2

  end subroutine quadratic_factors


  subroutine quadratic_factors_deriv(r, rj_lo, rj, rj_hi, dh_quadratic_dr)
    integer, intent(in)           ::  r, rj_lo, rj, rj_hi
    real, intent(out)             ::  dh_quadratic_dr(4)

    dh_quadratic_dr(1) = 4.0d0 * (-2.0d0*r + rj + rj_lo) / (rj - rj_lo)**2
    dh_quadratic_dr(2) = 0.0d0
    dh_quadratic_dr(3) = (4.0d0*r - rj - 3.0d0*rj_lo) / (rj - rj_lo)**2
    dh_quadratic_dr(4) = (4.0d0*r - rj - 3.0d0*rj_hi) / (rj_hi - rj)**2

  end subroutine quadratic_factors_deriv


  subroutine cubic_factors(r, rj_lo, rj, rj_hi, h_cubic)
    integer, intent(in)           :: r, rj_lo, rj, rj_hi
    real, intent(out)             :: h_cubic(4)

    h_cubic(1) =  3.0d0 * ( (r - rj_lo) / (rj - rj_lo) )**2 &
                 -2.0d0 * ( (r - rj_lo) / (rj - rj_lo) )**3
    h_cubic(2) =  3.0d0 * ( (rj_hi - r) / (rj_hi - rj) )**2 &
                 -2.0d0 * ( (rj_hi - r) / (rj_hi - rj) )**3
    h_cubic(3) = (r - rj) * ( (r - rj_lo) / (rj - rj_lo) )**2
    h_cubic(4) = (r - rj) * ( (r - rj_hi) / (rj_hi - rj) )**2

  end subroutine cubic_factors


  subroutine cubic_factors_deriv(r, rj_lo, rj, rj_hi, dh_cubic_dr)
    integer, intent(in)           :: r, rj_lo, rj, rj_hi
    real, intent(out)             :: dh_cubic_dr(4)

    dh_cubic_dr(1) =  6.0d0 * (r - rj_lo) / (rj - rj_lo)**2 &
                     -6.0d0 * (r - rj_lo)**2 / (rj - rj_lo)**3
    dh_cubic_dr(2) = -6.0d0 * (rj_hi - r) / (rj_hi - rj)**2 &
                     +6.0d0 * (rj_hi - r)**2 / (rj_hi - rj)**3
    dh_cubic_dr(3) = ( 2.0d0*(r - rj) * (r - rj_lo) + (r - rj_lo)**2 ) &
                     / (rj - rj_lo)**2
    dh_cubic_dr(4) = ( 2.0d0*(r - rj) * (r - rj_hi) + (r - rj_hi)**2 ) &
                     / (rj_hi - rj)**2

  end subroutine cubic_factors_deriv





end module mod_spline_functions
