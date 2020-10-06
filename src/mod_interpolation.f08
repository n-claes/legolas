module mod_interpolation
  use mod_global_variables, only: dp, ncool
  implicit none

  private

  public :: interpolate_table

contains


  !> Interpolates a given set of tables (x, y(x)) into a smooth curve.
  !! Assumes that x_table is an array with a monotone increase in values.
  !! Interpolation is done using <tt>ncool</tt> points, in general a second
  !! order polynomial approximation is used except near sharp jumps.
  !! @warning Throws an error if <tt>x_table</tt> is not monotone. @endwarning
  subroutine interpolate_table(n_table, x_table, y_table, x_interp, y_interp)
    use mod_logging, only: log_message

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

    ! check if x_table is monotonely increasing
    do i = 1, size(x_table)
      if (x_table(i + 1) < x_table(i)) then
        call log_message( &
          "interpolation: x-values are not monotonely increasing!", level="error" &
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

end module mod_interpolation
