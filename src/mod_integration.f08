! =============================================================================
!> Module responsible for integration of differential equations, useful
!! when setting equilibria or integrating the equilibrium equation.
!! Contains subroutines to numerically solve the following systems
!! of differential equations:
!! $$ y'(x) = A(x)y(x) + B(x) $$
!! These are solved using a fifth-order Runge-Kutta method.
module mod_integration
  use mod_global_variables, only: dp, dp_LIMIT
  implicit none

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
    xvalues, &
    axvalues, &
    bxvalues, &
    nbpoints, &
    yinit, &
    yvalues, &
    adaptive, &
    epsilon, &
    max_step_change, &
    min_dh_size, &
    max_dh_size, &
    max_iter_per_step, &
    new_xvalues &
  )
    use mod_interpolation, only: interpolate_table
    use mod_logging, only: log_message, str

    !> array of x-values
    real(dp), intent(in)  :: xvalues(:)
    !> term \(A(x)\)
    real(dp), intent(in)  :: axvalues(:)
    !> term \(B(x)\)
    real(dp), intent(in)  :: bxvalues(:)
    !> desired resolution of coefficient arrays <tt>axvalues</tt> and <tt>bxvalues</tt>
    integer, intent(in)   :: nbpoints
    !> initial value for \(y\) at left edge (<tt>xvalues(1)</tt>)
    real(dp), intent(in)  :: yinit
    !> array of y-values, size will be <tt>nbpoints</tt> if <tt>adaptive=.false.</tt>,
    !! size will vary if <tt>adaptive=.true.</tt>
    real(dp), intent(out), allocatable :: yvalues(:)
    !> if <tt>.true.</tt> an adaptive stepsize will be used
    logical, intent(in), optional   :: adaptive
    !> maximum truncation error, determines adaptive step
    real(dp), intent(in), optional  :: epsilon
    !> maximum factor by which dh can change if step is adaptive
    real(dp), intent(in), optional  :: max_step_change
    !> limit on the minimal stepsize
    real(dp), intent(in), optional  :: min_dh_size
    !> limit on the maximal stepsize
    real(dp), intent(in), optional  :: max_dh_size
    !> limit on maximum number of iterations for a single step
    integer, intent(in), optional   :: max_iter_per_step
    !> array of x-values corresponding to <tt>yvalues</tt>
    real(dp), intent(out), allocatable, optional  :: new_xvalues(:)

    integer   :: i, re_iter, maxiter
    real(dp)  :: ax_interp(nbpoints), bx_interp(nbpoints), x_interp(nbpoints)
    real(dp)  :: xi, xend, yi, dh, dh_new, dh_min, dh_max
    real(dp)  :: ysolrk4, ysolrk5, err, tol, max_dh_fact
    logical   :: use_adaptive_stepping
    real(dp), allocatable :: xtemp(:), ytemp(:)

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

    ! determine starting stepsize
    dh = x_interp(2) - x_interp(1)
    ! set max dh change
    max_dh_fact = 1.5d0
    if (present(max_step_change)) then
      max_dh_fact = max_step_change
    end if
    ! check adaptive stepping tolerance
    tol = 1.0d-4
    if (present(epsilon)) then
      tol = epsilon
    end if
    ! check min/max stepsizes
    dh_min = 1.0d-10
    dh_max = dh
    maxiter = 50
    if (present(min_dh_size)) then
      dh_min = min_dh_size
    end if
    if (present(max_dh_size)) then
      dh_max = max_dh_size
    end if
    if (present(max_iter_per_step)) then
      maxiter = max_iter_per_step
    end if
    ! set adaptive stepping
    use_adaptive_stepping = .false.
    if (present(adaptive)) then
      use_adaptive_stepping = adaptive
      if (use_adaptive_stepping) then
        call log_message( &
          "integrating using adaptive stepping, max dh change = " &
            // str(max_dh_fact, fmt="f5.1"), &
          level="info" &
        )
      else
        call log_message( &
          "integrating using regular stepping, dh = " // str(dh), level="info" &
        )
      end if
    end if

    ! allocate temporary arrays, this should be more than sufficient
    if (use_adaptive_stepping) then
      allocate(xtemp(10 * nbpoints), ytemp(10 * nbpoints))
    else
      allocate(xtemp(nbpoints), ytemp(nbpoints))
    end if

    ! do integration
    i = 1
    xi = x_interp(1)
    xend = x_interp(nbpoints)
    xtemp(1) = xi
    ytemp(1) = yinit
    re_iter = 0
    do while (xi < xend)
      dh = min(dh, xend - xi)
      yi = ytemp(i)

      ! sanity check, for non-adaptive stepping we sometimes we have xend - xi
      ! equal to about 1e-13 at the last iteration i == nbpoints, and the above part
      ! selects that as dh near the end. This results in an "additional" iteration
      ! with a bogus dh due to numerical round-off, and a final array of length
      ! nbpoints + 1. Hence, if we're not using adaptive stepping and we encounter
      ! something like this, modify dh so that it nicely reaches xend next iteration.
      if ( &
        .not. use_adaptive_stepping &
        .and. i == nbpoints - 1 &
      ) then
        dh = xend - xi
      end if

      call rk45(xi, dh, x_interp, ax_interp, bx_interp, yi, ysolrk4, ysolrk5)
      err = abs(ysolrk5 - ysolrk4) / dh

      if (use_adaptive_stepping) then
        ! calculate truncation error and new stepsize
        dh_new = 0.9d0 * dh * (tol / err) ** 0.2d0
        ! make sure stepsize does not change too much
        if (dh_new / dh > max_dh_fact) then
          dh = dh * max_dh_fact     ! stepsize was increased too much
        else if (dh / dh_new > max_dh_fact) then
          dh = dh / max_dh_fact     ! stepsize was decreased too much
        else
          dh = dh_new               ! stepsize is fine, replace
        end if
        ! make sure stepsize is within bounds
        if (dh < dh_min) then
          dh = dh_min
        else if (dh > dh_max) then
          dh = dh_max
        end if
        ! if error is not within tolerance, retake this step with new dh
        ! continue if iterations is higher than allowed per step
        if (err > tol .and. re_iter < maxiter) then
          re_iter = re_iter + 1
          cycle
        end if
        re_iter = 0
        ! also handle tiny dh for adaptive stepping
        if ((xend - xi - dh) < dh_min) then
          dh = xend - xi
        end if
      end if

      ! we are either within tolerance or without adaptive stepping, take next step
      i = i + 1
      xi = xi + dh
      xtemp(i) = xi
      ytemp(i) = ysolrk5

      if (mod(i, 10000) == 0) then
        call log_message( &
          "steps: " // str(i) // " | xi: " // str(xi) // " | dh: " &
          // str(dh, fmt="e20.5"), &
          level="info", &
          use_prefix=.false. &
        )
      end if
    end do

    ! allocate final arrays
    allocate(yvalues(i))
    yvalues = ytemp(:i)
    deallocate(ytemp)
    if (present(new_xvalues)) then
      allocate(new_xvalues(i))
      new_xvalues = xtemp(:i)
      deallocate(xtemp)
    end if
  end subroutine integrate_ode_rk


  !> Calculates the Runge-Kutta coefficients and calculates the fourth and fifth
  !! order solutions for step i+1 based on the values at step i.
  subroutine rk45(xi, dh, xarray, axarray, bxarray, yi, ysolrk4, ysolrk5)
    use mod_interpolation, only: lookup_table_value
    !> current x value in iteration
    real(dp), intent(in)  :: xi
    !> current step size
    real(dp), intent(in)  :: dh
    !> array of x values
    real(dp), intent(in)  :: xarray(:)
    !> array of A(x) values
    real(dp), intent(in)  :: axarray(:)
    !> array of B(x) values
    real(dp), intent(in)  :: bxarray(:)
    !> current y value in iteration
    real(dp), intent(in)  :: yi
    !> fourth-order solution
    real(dp), intent(out) :: ysolrk4
    !> fifth-order solution
    real(dp), intent(out) :: ysolrk5

    real(dp)  :: axi, bxi, yvalrk
    real(dp)  :: rkf1, rkf2, rkf3, rkf4, rkf5, rkf6

    ! first step
    axi = lookup_table_value(xi, xarray, axarray)
    bxi = lookup_table_value(xi, xarray, bxarray)
    yvalrk = yi
    rkf1 = dh * (axi * yvalrk + bxi)
    ! second step
    axi = lookup_table_value(xi + 0.25d0 * dh, xarray, axarray)
    bxi = lookup_table_value(xi + 0.25d0 * dh, xarray, bxarray)
    yvalrk = yi + 0.25d0 * rkf1
    rkf2 = dh * (axi * yvalrk + bxi)
    ! third step
    axi = lookup_table_value(xi + 3.0d0 * dh / 8.0d0, xarray, axarray)
    bxi = lookup_table_value(xi + 3.0d0 * dh / 8.0d0, xarray, bxarray)
    yvalrk = yi + (3.0d0 * rkf1 + 9.0d0 * rkf2) / 32.0d0
    rkf3 = dh * (axi * yvalrk + bxi)
    ! fourth step
    axi = lookup_table_value(xi + 12.0d0 * dh / 13.0d0, xarray, axarray)
    bxi = lookup_table_value(xi + 12.0d0 * dh / 13.0d0, xarray, bxarray)
    yvalrk = yi + (1932.0d0 * rkf1 - 7200.0d0 * rkf2 + 7296.0d0 * rkf3) / 2197.0d0
    rkf4 = dh * (axi * yvalrk + bxi)
    ! fifth step
    axi = lookup_table_value(xi + dh, xarray, axarray)
    bxi = lookup_table_value(xi + dh, xarray, bxarray)
    yvalrk = ( &
      yi &
      + (439.0d0 / 216.0d0) * rkf1 &
      - 8.0d0 * rkf2 &
      + (3680.0d0 / 513.0d0) * rkf3 &
      - (845.0d0 / 4104.0d0) * rkf4 &
    )
    rkf5 = dh * (axi * yvalrk + bxi)
    ! sixth step
    axi = lookup_table_value(xi + 0.5d0 * dh, xarray, axarray)
    bxi = lookup_table_value(xi + 0.5d0 * dh, xarray, bxarray)
    yvalrk = ( &
    yi &
      - (8.0d0 / 27.0d0) * rkf1 &
      + 2.0d0 * rkf2 &
      - (3544.0d0 / 2565.0d0) * rkf3 &
      + (1859.0d0 / 4104.0d0) * rkf4 &
      - (11.0d0 / 40.0d0) * rkf5 &
    )
    rkf6 = dh * (axi * yvalrk + bxi)

    ! fourth order solution
    ysolrk4 = ( &
      yi &
      + (25.0d0 / 216.0d0) * rkf1 &
      + (1408.0d0 / 2565.0d0) * rkf3 &
      + (2197.0d0 / 4104.0d0) * rkf4 &
      - (1.0d0 / 5.0d0) * rkf5 &
    )
    ! fifth order solution
    ysolrk5 = ( &
      yi &
      + (16.0d0 / 135.0d0) * rkf1 &
      + (6656.0d0 / 12825.0d0) * rkf3 &
      + (28561.0d0 / 56430.0d0) * rkf4 &
      - (9.0d0 / 50.0d0) * rkf5 &
      + (2.0d0 / 55.0d0) * rkf6 &
    )
  end subroutine rk45


  !> Checks if an array needs resampling.
  function needs_resampling(base, target, array_name) result(resample)
    use mod_logging, only: log_message, str

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
      call log_message( &
        "ode integrator: resampling " // trim(adjustl(array_name)) // " (" &
        // str(base) // " -> " // str(target) // " points)", &
        level="info" &
      )
      resample = .true.
    end if
  end function needs_resampling
end module mod_integration
