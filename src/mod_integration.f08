! =============================================================================
!> Module responsible for integration of differential equations, useful
!! when setting equilibria or integrating the equilibrium equation.
!! Contains subroutines to numerically solve the following systems
!! of differential equations:
!! $$ y'(x) = A(x)y(x) + B(x) $$
!! These are solved using a fifth-order Runge-Kutta method.
module mod_integration
  use mod_global_variables, only: dp
  use mod_logging, only: logger, str, exp_fmt
  implicit none

  private

  interface
    real(dp) function func(x)
      use mod_global_variables, only: dp
      real(dp), intent(in) :: x
    end function func
  end interface

  public :: integrate_ode_rk45

contains

  !> Numerically integrates the differential equation
  !! $$ y'(x) = A(x)y(x) + B(x) $$ using a fifth-order Runge-Kutta method.
  !! The functions A(x) and B(x) are passed as arguments and should be conform to the
  !! given interface, that is, these should be `real(dp)` functions which take a single
  !! `real(dp), intent(in)` argument.
  !! The integration is performed on the interval [x0, x1] with a stepsize of
  !! `dh = (x1 - x0) / (nbpoints - 1)`.
  subroutine integrate_ode_rk45( &
    x0, x1, ax_func, bx_func, nbpoints, yinit, yvalues, xvalues &
  )
    !> start of x-interval
    real(dp), intent(in) :: x0
    !> end of x-interval
    real(dp), intent(in) :: x1
    !> function to calculate A(x)
    procedure(func) :: ax_func
    !> function to calculate B(x)
    procedure(func) :: bx_func
    !> number of points, determines stepsize
    integer, intent(in) :: nbpoints
    !> initial value of y
    real(dp), intent(in) :: yinit
    !> integrated values of y
    real(dp), intent(out), allocatable :: yvalues(:)
    !> x-values corresponding to integrated y-values
    real(dp), intent(out), allocatable, optional :: xvalues(:)

    real(dp) :: xi, yi, dh
    real(dp) :: ysolrk4, ysolrk5
    real(dp), allocatable :: xivalues(:)
    integer :: i

    dh = (x1 - x0) / (nbpoints - 1)
    call logger%info( &
      "integrating from " // str(x0) // " to " // str(x1) &
      // " with dh = " // str(dh, fmt=exp_fmt) &
    )
    xi = x0
    yi = yinit
    allocate(yvalues(nbpoints), xivalues(nbpoints))
    xivalues(1) = xi
    yvalues(1) = yi
    do i = 2, nbpoints - 1
      call rk45(xi, dh, ax_func, bx_func, yi, ysolrk4, ysolrk5)
      yi = ysolrk5
      xi = xi + dh
      yvalues(i) = yi
      xivalues(i) = xi
    end do
    ! final step, do outside loop to ensure dh ends up at x1 (rounding errors)
    dh = x1 - xi
    call rk45(xi, dh, ax_func, bx_func, yi, ysolrk4, ysolrk5)
    yvalues(nbpoints) = ysolrk5
    xivalues(nbpoints) = x1
    if (present(xvalues)) xvalues = xivalues
    deallocate(xivalues)
  end subroutine integrate_ode_rk45


  !> Calculates the Runge-Kutta coefficients and calculates the fourth and fifth
  !! order solutions for step i+1 based on the values at step i.
  subroutine rk45(xi, dh, ax_func, bx_func, yi, ysolrk4, ysolrk5)
    !> current x value
    real(dp), intent(in) :: xi
    !> current step size
    real(dp), intent(in) :: dh
    !> function to calculate A(x)
    procedure(func) :: ax_func
    !> function to calculate B(x)
    procedure(func) :: bx_func
    !> current y value
    real(dp), intent(in) :: yi
    !> fourth order solution
    real(dp), intent(out) :: ysolrk4
    !> fifth order solution
    real(dp), intent(out) :: ysolrk5

    real(dp) :: rkf1, rkf2, rkf3, rkf4, rkf5, rkf6
    real(dp) :: xvalrk, yvalrk

    ! first step
    rkf1 = dh * (ax_func(xi) * yi + bx_func(xi))
    ! second step
    xvalrk = xi + 0.25_dp * dh
    yvalrk = yi + 0.25_dp * rkf1
    rkf2 = dh * (ax_func(xvalrk) * yvalrk + bx_func(xvalrk))
    ! third step
    xvalrk = xi + 3.0_dp * dh / 8.0_dp
    yvalrk = yi + (3.0_dp * rkf1 + 9.0_dp * rkf2) / 32.0_dp
    rkf3 = dh * (ax_func(xvalrk) * yvalrk + bx_func(xvalrk))
    ! fourth step
    xvalrk = xi + 12.0_dp * dh / 13.0_dp
    yvalrk = yi + (1932.0_dp * rkf1 - 7200.0_dp * rkf2 + 7296.0_dp * rkf3) / 2197.0_dp
    rkf4 = dh * (ax_func(xvalrk) * yvalrk + bx_func(xvalrk))
    ! fifth step
    xvalrk = xi + dh
    yvalrk = yi + ( &
      439.0_dp * rkf1 / 216.0_dp &
      - 8.0_dp * rkf2 &
      + 3680.0_dp * rkf3 / 513.0_dp &
      - 845.0_dp * rkf4 / 4104.0_dp &
    )
    rkf5 = dh * (ax_func(xvalrk) * yvalrk + bx_func(xvalrk))
    ! sixth step
    xvalrk = xi + 0.5_dp * dh
    yvalrk = yi + ( &
      - 8.0_dp * rkf1 / 27.0_dp &
      + 2.0_dp * rkf2 &
      - 3544.0_dp * rkf3 / 2565.0_dp &
      + 1859.0_dp * rkf4 / 4104.0_dp &
      - 11.0_dp * rkf5 / 40.0_dp &
    )
    rkf6 = dh * (ax_func(xvalrk) * yvalrk + bx_func(xvalrk))

    ! fourth order solution
    ysolrk4 = yi + ( &
      25.0_dp * rkf1 / 216.0_dp &
      + 1408.0_dp * rkf3 / 2565.0_dp &
      + 2197.0_dp * rkf4 / 4104.0_dp &
      - rkf5 / 5.0_dp &
    )
    ! fifth order solution
    ysolrk5 = yi + ( &
      16.0_dp * rkf1 / 135.0_dp &
      + 6656.0_dp * rkf3 / 12825.0_dp &
      + 28561.0_dp * rkf4 / 56430.0_dp &
      - 9.0_dp * rkf5 / 50.0_dp &
      + 2.0_dp * rkf6 / 55.0_dp &
    )
  end subroutine rk45


end module mod_integration
