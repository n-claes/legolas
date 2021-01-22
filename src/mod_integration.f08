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
  !! method. The argument <tt>nbpoints</tt> determines the stepsize through
  !! $$ dh = \frac{xvalues(N) - xvalues(1)}{nbpoints} $$
  !! If the arrays \(A(x), B(x)\) are not of size <tt>nbpoints</tt>, then these are
  !! interpolated to that resolution. The differential equation is then integrated,
  !! the solution will also be of size <tt>nbpoints</tt> and can be downsampled using
  !! the appropriate subroutine.
  !! If desired the optional argument <tt>dyvalues</tt> can be provided,
  !! which contains the (numerical) derivative of y.
  subroutine integrate_ode_rk( &
    xvalues, axvalues, bxvalues, nbpoints, yinit, yvalues, dyvalues &
  )
    use mod_interpolation, &
      only: interpolate_table, lookup_table_value, get_numerical_derivative
    use mod_logging, only: log_message, char_log

    !> array of x-values
    real(dp), intent(in)  :: xvalues(:)
    !> term \(A(x)\)
    real(dp), intent(in)  :: axvalues(:)
    !> term \(B(x)\)
    real(dp), intent(in)  :: bxvalues(:)
    !> number of points to do the integration
    integer, intent(in)   :: nbpoints
    !> initial value for \(y\) at left edge (<tt>xvalues(1)</tt>)
    real(dp), intent(in)  :: yinit
    !> array of y-values, of length <tt>nbpoints</tt>
    real(dp), intent(out) :: yvalues(nbpoints)
    !> array of dy/dx-values
    real(dp), intent(out), optional :: dyvalues(nbpoints)

    integer   :: i
    real(dp)  :: ax_interp(nbpoints), bx_interp(nbpoints), x_interp(nbpoints)
    real(dp)  :: xi, yi, axi, bxi, yvalrk

    ! check if A(x) needs to be resampled
    if (needs_resampling(size(axvalues), nbpoints, "A(x)")) then
      call interpolate_table(nbpoints, xvalues, axvalues, x_interp, ax_interp)
    else
      ax_interp = axvalues
      x_interp = xvalues
    end if
    ! check if B(x) needs to be resampled
    if (needs_resampling(size(bxvalues), nbpoints, "B(x)")) then
      call interpolate_table(nbpoints, xvalues, bxvalues, x_interp, bx_interp)
    else
      bx_interp = bxvalues
      x_interp = xvalues
    end if

    ! do actual integration, first fill initial value
    yvalues(1) = yinit
    do i = 1, nbpoints - 1
      xi = x_interp(i)
      yi = yvalues(i)
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
      yvalrk = ( &
        yi &
        + (439.0d0 / 216.0d0) * rkf1 &
        - 8.0d0 * rkf2 &
        + (3680.0d0 / 513.0d0) * rkf3 &
        - (845.0d0 / 4104.0d0) * rkf4 &
      )
      rkf5 = dh * (axi * yvalrk + bxi)

      ! sixth step
      axi = lookup_table_value(xi + 0.5d0 * dh, x_interp, ax_interp)
      bxi = lookup_table_value(xi + 0.5d0 * dh, x_interp, bx_interp)
      yvalrk = ( &
        yi &
        - (8.0d0 / 27.0d0) * rkf1 &
        + 2.0d0 * rkf2 &
        - (3544.0d0 / 2565.0d0) * rkf3 &
        + (1859.0d0 / 4104.0d0) * rkf4 &
        - (11.0d0 / 40.0d0) * rkf5 &
      )
      rkf6 = dh * (axi * yvalrk + bxi)

      yvalues(i + 1) = ( &
        yi + (16.0d0 / 135.0d0) * rkf1 &
        + (6656.0d0 / 12825.0d0) * rkf3 &
        + (28561.0d0 / 56430.0d0) * rkf4 &
        - (9.0d0 / 50.0d0) * rkf5 &
        + (2.0d0 / 55.0d0) * rkf6 &
      )

      if (mod(i - 1, int(nbpoints / 10)) == 0) then
        write(char_log, '(f20.2)') i * 100.0d0 / nbpoints
        call log_message( &
          "integration progress: " // trim(adjustl(char_log)) // "%", &
          level="info" &
        )
      end if
    end do
    call log_message("integration progress: 100%", level="info")

    ! fill optional return array
    if (present(dyvalues)) then
      call get_numerical_derivative(x_interp, yvalues, dyvalues)
    end if
  end subroutine integrate_ode_rk


  !> Checks if an array needs resampling.
  function needs_resampling(base, target, array_name) result(resample)
    use mod_logging, only: log_message, char_log, char_log2, int_fmt
    !> size of base array
    integer, intent(in) :: base
    !> size of target array
    integer, intent(in) :: target
    !> name of array (printed to console)
    character(len=*), intent(in)  :: array_name
    !> is <tt>.true.</tt> if arrays differ in size
    logical :: resample

    if (target < base) then
      call log_message( &
        "resampling: target #points is less than size of input arrays!", &
        level="warning" &
      )
    end if

    if (base == target) then
      resample = .false.
    else
      write(char_log, int_fmt) base
      write(char_log2, int_fmt) target
      call log_message( &
        "ode integrator: resampling " // trim(adjustl(array_name)) // " (" &
        // trim(adjustl(char_log)) &
        // " -> " // trim(adjustl(char_log2)) // " points)", &
        level="info" &
      )
      resample = .true.
    end if
  end function needs_resampling
end module mod_integration
