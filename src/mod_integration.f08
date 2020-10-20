! =============================================================================
!> Module responsible for integration of differential equations, useful
!! when setting equilibria or integrating the equilibrium equation.
!! Contains subroutines to numerically solve the following systems
!! of differential equations:
!! $$ y'(x) = A(x)y(x) + B(x) $$
!! These are solved using a fifth-order Runge-Kutta method.
module mod_integration
  use mod_global_variables, only: dp
  implicit none

  !> step size for discretisation
  real(dp)  :: dh
  !> first Runge-Kutta factor
  real(dp)  :: rkf1
  !> second Runge-Kutta factor
  real(dp)  :: rkf2
  !> third Runge-Kutta factor
  real(dp)  :: rkf3
  !> fourth Runge-Kutta factor
  real(dp)  :: rkf4
  !> fifth Runge-Kutta factor
  real(dp)  :: rkf5
  !> sixth Runge-Kutta factor
  real(dp)  :: rkf6

  private

  public :: integrate_ode_rk

contains


  !> Integrates a first order differential equation of the form
  !! $$ y'(x) = A(x)y(x) + B(x) $$ using a fifth-order Runge-Kutta
  !! method. Both \(A(x), B(x)\) and <tt>xvalues</tt> are interpolated at high
  !! (<tt>nbpoints</tt>) resolution. The equation is then integrated at this
  !! resolution, after which the high-res solution is downsampled to be of the
  !! same size as the input arrays <tt>xvalues</tt>.
  subroutine integrate_ode_rk(xvalues, axvalues, bxvalues, yvalues, yinit, nbpoints)
    use mod_interpolation, only: interpolate_table, lookup_table_value
    use mod_logging, only: log_message


    !> array of x-values
    real(dp), intent(in)  :: xvalues(:)
    !> term \(A(x)\)
    real(dp), intent(in)  :: axvalues(:)
    !> term \(B(x)\)
    real(dp), intent(in)  :: bxvalues(:)
    !> array of (integrated) y-values, same size as <tt>xvalues</tt>
    real(dp), intent(out) :: yvalues(size(xvalues))
    !> initial value for \(y\) at left edge (<tt>xvalues(1)</tt>)
    real(dp), intent(in)  :: yinit
    !> number of points to do the integration
    integer, intent(in)   :: nbpoints

    integer   :: i
    real(dp)  :: ax_interp(nbpoints), bx_interp(nbpoints), x_interp(nbpoints)
    real(dp)  :: yvalues_hr(nbpoints)
    real(dp)  :: xi, yi, axi, bxi, yvalrk

    if (nbpoints < size(xvalues)) then
      call log_message( &
        "integrate_ode_rk: nbpoints is smaller than size of input arrays!", &
        level="warning" &
      )
    end if

    ! fill initial value
    yvalues_hr(1) = yinit

    call interpolate_table(nbpoints, xvalues, axvalues, x_interp, ax_interp)
    call interpolate_table(nbpoints, xvalues, bxvalues, x_interp, bx_interp)

    do i = 1, nbpoints - 1
      xi = x_interp(i)
      yi = yvalues_hr(i)
      dh = x_interp(i + 1) - x_interp(i)

      ! first step
      axi = lookup_table_value(xi, x_interp, ax_interp)
      bxi = lookup_table_value(xi, x_interp, bx_interp)
      yvalrk = yi
      rkf1 = dh * (axi * yvalrk + bxi)

      ! second step
      axi = lookup_table_value(xi + 0.25d0 * dh, x_interp, ax_interp)
      bxi = lookup_table_value(xi + 0.25d0 * dh, x_interp, bx_interp)
      yvalrk = yi + 0.25d0 * rkf1
      rkf2 = dh * (axi * yvalrk + bxi)

      ! third step
      axi = lookup_table_value(xi + 3.0d0 * dh / 8.0d0, x_interp, ax_interp)
      bxi = lookup_table_value(xi + 3.0d0 * dh / 8.0d0, x_interp, bx_interp)
      yvalrk = yi + (3.0d0 * rkf1 + 9.0d0 * rkf2) / 32.0d0
      rkf3 = dh * (axi * yvalrk + bxi)

      ! fourth step
      axi = lookup_table_value(xi + 12.0d0 * dh / 13.0d0, x_interp, ax_interp)
      bxi = lookup_table_value(xi + 12.0d0 * dh / 13.0d0, x_interp, bx_interp)
      yvalrk = yi + (1932.0d0 * rkf1 - 7200.0d0 * rkf2 + 7296.0d0 * rkf3) / 2197.0d0
      rkf4 = dh * (axi * yvalrk + bxi)

      ! fifth step
      axi = lookup_table_value(xi + dh, x_interp, ax_interp)
      bxi = lookup_table_value(xi + dh, x_interp, bx_interp)
      yvalrk = yi + (439.0d0 / 216.0d0) * rkf1 - 8.0d0 * rkf2 &
               + (3680.0d0 / 513.0d0) * rkf3 - (845.0d0 / 4104.0d0) * rkf4
      rkf5 = dh * (axi * yvalrk + bxi)

      ! sixth step
      axi = lookup_table_value(xi + 0.5d0 * dh, x_interp, ax_interp)
      bxi = lookup_table_value(xi + 0.5d0 * dh, x_interp, bx_interp)
      yvalrk = yi - (8.0d0 / 27.0d0) * rkf1 + 2.0d0 * rkf2 &
               - (3544.0d0 / 2565.0d0) * rkf3 + (1859.0d0 / 4104.0d0) * rkf4 &
               - (11.0d0 / 40.0d0) * rkf5
      rkf6 = dh * (axi * yvalrk + bxi)

      yvalues_hr(i + 1) = yi + (16.0d0 / 135.0d0) * rkf1 &
                          + (6656.0d0 / 12825.0d0) * rkf3 &
                          + (28561.0d0 / 56430.0d0) * rkf4 &
                          - (9.0d0 / 50.0d0) * rkf5 &
                          + (2.0d0 / 55.0d0) * rkf6
    end do

    ! fill return array
    do i = 1, size(yvalues)
      yvalues(i) = lookup_table_value(xvalues(i), x_interp, yvalues_hr)
    end do
  end subroutine integrate_ode_rk
end module mod_integration