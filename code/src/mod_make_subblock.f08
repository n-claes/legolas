module mod_make_subblock
  use mod_global_variables
  implicit none


contains

    subroutine subblock(quadblock, factors, len_factors, spline1, spline2)

      complex(dp), intent(in) :: quadblock(dim_quadblock, dim_quadblock)
      complex(dp), intent(in) :: factors(5)
      integer, intent(in)     :: len_factors
      real(dp), intent(in)    :: spline1(4), spline2(4)
      return
    end subroutine subblock


end module mod_make_subblock
