module mod_make_subblock
  use mod_global_variables
  implicit none


contains

    subroutine subblock(quadblock, factors, positions, &
                        curr_weight, spline1, spline2)

      complex(dp), intent(in)     :: quadblock(dim_quadblock, dim_quadblock)
      complex(dp), intent(inout)  :: factors(:)
      integer, intent(in)         :: positions(:, :)
      real(dp), intent(in)        :: curr_weight, spline1(4), spline2(4)

      integer                     :: i, len_factors

      len_factors = size(factors)

      do i = 1, len_factors
        factors(i) = factors(i) * curr_weight
      end do

    end subroutine subblock


end module mod_make_subblock
