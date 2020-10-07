! =============================================================================
!> Module responsible for table interpolations and array lookups.
!! Contains subroutines for table interpolations, numerical derivatives
!! of arrays and lookup functions.
!! Subroutines are loosely based on routines implemented in the
!! [MPI-AMRVAC](amrvac.org) code.
module mod_interpolation
  use mod_global_variables, only: dp, ncool
  use mod_logging, only: log_message
  implicit none

  private

  public :: interpolate_table
  public :: get_numerical_derivative
  public :: lookup_table_value

contains


  !> Interpolates a given set of tables (x, y(x)) into a smooth curve.
  !! Assumes that x_table is an array with a monotone increase in values.
  !! Interpolation is done using <tt>ncool</tt> points, in general a second
  !! order polynomial approximation is used except near sharp jumps.
  !! @warning Throws an error if <tt>x_table</tt> is not monotone. @endwarning
  subroutine interpolate_table(n_table, x_table, y_table, x_interp, y_interp)
    !> number of points in the tables
    integer, intent(in)   :: n_table
    !> x-values in the table
    real(dp), intent(in)  :: x_table(:)
    !> y-values in the table
    real(dp), intent(in)  :: y_table(:)
    !> interpolated x-values
    real(dp), intent(out) :: x_interp(ncool)
    !> interpolated y-values
    real(dp), intent(out) :: y_interp(ncool)

    integer   :: i, j
    real(dp)  :: fact1, fact2, fact3
    real(dp)  :: dx, dy1, dy2
    logical   :: jump

    ! check if x_table is a monotonically increasing array
    do i = 1, size(x_table) - 1
      if (x_table(i + 1) < x_table(i)) then
        call log_message( &
          "interpolation: x-values are not monotonically increasing!", level="error" &
        )
      end if
    end do

    ! outer edges remain the same
    x_interp(1) = x_table(1)
    x_interp(ncool) = x_table(n_table)
    y_interp(1) = y_table(1)
    y_interp(ncool) = y_table(n_table)

    dx = (x_table(n_table) - x_table(1)) / (ncool - 1)
    do i = 2, ncool
      x_interp(i) = x_interp(i - 1) + dx
      do j = 1, n_table - 1
        jump = .false.
        if (x_interp(i) < x_table(j + 1)) then
          if (j < n_table - 1) then
            ! check if we have a sharp jump
            dy1 = y_table(j + 1) - y_table(j)
            dy2 = y_table(j + 2) - y_table(j + 1)
            jump = ( max(dabs(dy1), dabs(dy2)) > 2 * min(dabs(dy1), dabs(dy2)) )
          end if
          ! no interpolation at outer edge or near sharp jumps
          if ((j == n_table - 1) .or. jump) then
            fact1 = (x_interp(i) - x_table(j + 1)) / (x_table(j) - x_table(j + 1))
            fact2 = (x_interp(i) - x_table(j)) / (x_table(j + 1) - x_table(j))
            y_interp(i) = fact1 * y_table(j) + fact2 * y_table(j + 1)
            exit
          end if

          fact1 = ((x_interp(i) - x_table(j + 1)) * (x_interp(i) - x_table(j + 2))) &
                / ((x_table(j) - x_table(j + 1)) * (x_table(j) - x_table(j + 2)))
          fact2 = ((x_interp(i) - x_table(j)) / (x_interp(i) - x_table(j + 2))) &
                / ((x_table(j + 1) - x_table(j)) * (x_table(j + 1) - x_table(j + 2)))
          fact3 = ((x_interp(i) - x_table(j)) * (x_interp(i) - x_table(j + 1))) &
                / ((x_table(j + 2) - x_table(j)) * (x_table(j + 2) - x_table(j + 1)))
          y_interp(i) = fact1 * y_table(j) + fact2 * y_table(j + 1) &
                      + fact3 * y_table(j + 2)
          exit
        end if
      end do
    end do
  end subroutine interpolate_table


  !> Calculates the numerical derivative of a given array.
  !! The outer values in the array use their immediate neighbours (index 2 and
  !! index n-1) to calculate the derivative. A simple two-point method is used,
  !! where the derivative at index \(i\) is given by
  !! $$ dy_i = \frac{y_{i+1} - y_{i-1}}{x_{i+1} - x_{i-1}} $$
  !! Returns an array of the same size as <tt>y_values</tt>.
  !! @warning Throws an error if <tt>x_values</tt> and <tt>y_values</tt> differ
  !!          in size. @endwarning
  subroutine get_numerical_derivative(x_values, y_values, dy_values)
    !> x-values against which to differentiate
    real(dp), intent(in)  :: x_values(:)
    !> array of y-values, assuming \(y(x)\) relation
    real(dp), intent(in)  :: y_values(:)
    !> derivative of \(y\) with respect to \(x\), same size as input arrays
    real(dp), intent(out) :: dy_values(size(y_values))

    integer :: i, nvals

    ! x_values and y_values should be the same length
    if (size(x_values) /= size(y_values)) then
      call log_message( &
        "numerical derivative: x_values and y_values should have the same size!", &
        level="error" &
      )
    end if

    nvals = size(x_values)
    dy_values(1) = (y_values(2) - y_values(1)) / (x_values(2) - x_values(1))
    dy_values(nvals) = (y_values(nvals) - y_values(nvals - 1)) &
                     / (x_values(nvals) - x_values(nvals - 1))
    do i = 2, nvals - 1
      dy_values(i) = (y_values(i + 1) - y_values(i - 1)) &
                   / (x_values(i + 1) - x_values(i - 1))
    end do
  end subroutine get_numerical_derivative


  !> Function for fast table-lookup, returns the corresponding y-value
  !! in <tt>y_values</t>> based on a given based on a given \(x0\).
  !! Uses simple linear interpolation.
  !! @warning Throws an error if
  !!
  !! - <tt>x_values</tt> is not monotonically increasing.
  !! - x is outside of <tt>x_values</tt> range. @endwarning
  function lookup_table_value(x, x_values, y_values) result(y_found)
    !> value to look up
    real(dp), intent(in)  :: x
    !> array of x-values
    real(dp), intent(in)  :: x_values(:)
    !> array of y-values, assuming \(y(x)\) relation
    real(dp), intent(in)  :: y_values(:)

    integer   :: idx, nvals
    real(dp)  :: x0, x1, y0, y1
    !> interpolated y-value based on \(x0\)
    real(dp)  :: y_found

    nvals = size(x_values)
    ! check if x_values is a monotonically increasing array
    do idx = 1, nvals - 1
      if (x_values(idx + 1) < x_values(idx)) then
        call log_message( &
          "lookup_value: x_values are not monotonically increasing!", level="error" &
        )
      end if
    end do

    ! check if we are outside of the table
    if (x < x_values(1)) then
      call log_message("lookup_value: x outside x_values (too small)", level="error")
    else if (x > x_values(nvals)) then
      call log_message("lookup_value: x outside x_values (too large)", level="error")
    end if

    ! index of nearest value to x (dim=1 to return a scalar)
    idx = minloc(abs(x_values - x), dim=1)

    ! check edges
    if (idx == 1) then
      y_found = y_values(1)
      return
    else if (idx == nvals) then
      y_found = y_values(nvals)
      return
    end if

    ! check if we're on left of right side of nearest point
    if (x < x_values(idx)) then
      x0 = x_values(idx - 1)
      x1 = x_values(idx)
      y0 = y_values(idx - 1)
      y1 = y_values(idx)
    else
      x0 = x_values(idx)
      x1 = x_values(idx + 1)
      y0 = y_values(idx)
      y1 = y_values(idx + 1)
    end if

    ! do linear interpolation
    y_found = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
  end function lookup_table_value

end module mod_interpolation
