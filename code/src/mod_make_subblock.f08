module mod_make_subblock
  use mod_global_variables
  implicit none

  public


contains

    subroutine subblock(quadblock, factors, positions, &
                        curr_weight, spline1, spline2)

      complex(dp), intent(inout)  :: quadblock(dim_quadblock, dim_quadblock)
      complex(dp), intent(inout)  :: factors(:)
      integer, intent(in)         :: positions(:, :)
      real(dp), intent(in)        :: curr_weight, spline1(4), spline2(4)

      integer                     :: i, len_factors
      integer                     :: curr_position(2)
      integer                     :: idx(2)

      idx(:) = 0

      len_factors = size(factors)

      do i = 1, len_factors
        factors(i) = factors(i) * curr_weight
      end do

      do i = 1, len_factors
        curr_position = positions(i, :)

        ! Subblock top-left corner
        idx(:) = curr_position(:)
        quadblock(idx(1)  , idx(2)  ) = quadblock(idx(1)  , idx(2)  ) + &
                                        spline1(2) * factors(i) * spline2(2)
        quadblock(idx(1)  , idx(2)+1) = quadblock(idx(1)  , idx(2)+1) + &
                                        spline1(2) * factors(i) * spline2(4)
        quadblock(idx(1)+1, idx(2)  ) = quadblock(idx(1)+1, idx(2)  ) + &
                                        spline1(4) * factors(i) * spline2(2)
        quadblock(idx(1)+1, idx(2)+1) = quadblock(idx(1)+1, idx(2)+1) + &
                                        spline1(4) * factors(i) * spline2(4)

        ! Subblock top-right corner
        idx(:) = [curr_position(1), curr_position(2) + dim_subblock]
        quadblock(idx(1)  , idx(2)  ) = quadblock(idx(1)  , idx(2)  ) + &
                                        spline1(2) * factors(i) * spline2(1)
        quadblock(idx(1)  , idx(2)+1) = quadblock(idx(1)  , idx(2)+1) + &
                                        spline1(2) * factors(i) * spline2(3)
        quadblock(idx(1)+1, idx(2)  ) = quadblock(idx(1)+1, idx(2)  ) + &
                                        spline1(4) * factors(i) * spline2(1)
        quadblock(idx(1)+1, idx(2)+1) = quadblock(idx(1)+1, idx(2)+1) + &
                                        spline1(4) * factors(i) * spline2(3)

        ! Subblock bottom-left corner
        idx(:) = [curr_position(1) + dim_subblock, curr_position(2)]
        quadblock(idx(1)  , idx(2)  ) = quadblock(idx(1)  , idx(2)  ) + &
                                        spline1(1) * factors(i) * spline2(2)
        quadblock(idx(1)  , idx(2)+1) = quadblock(idx(1)  , idx(2)+1) + &
                                        spline1(1) * factors(i) * spline2(4)
        quadblock(idx(1)+1, idx(2)  ) = quadblock(idx(1)+1, idx(2)  ) + &
                                        spline1(3) * factors(i) * spline2(2)
        quadblock(idx(1)+1, idx(2)+1) = quadblock(idx(1)+1, idx(2)+1) + &
                                        spline1(3) * factors(i) * spline2(4)

        ! Subblock bottom-right corner
        idx(:) = curr_position(:) + dim_subblock
        quadblock(idx(1)  , idx(2)  ) = quadblock(idx(1)  , idx(2)  ) + &
                                        spline1(1) * factors(i) * spline2(1)
        quadblock(idx(1)  , idx(2)+1) = quadblock(idx(1)  , idx(2)+1) + &
                                        spline1(1) * factors(i) * spline2(3)
        quadblock(idx(1)+1, idx(2)  ) = quadblock(idx(1)+1, idx(2)  ) + &
                                        spline1(3) * factors(i) * spline2(1)
        quadblock(idx(1)+1, idx(2)+1) = quadblock(idx(1)+1, idx(2)+1) + &
                                        spline1(3) * factors(i) * spline2(3)
      end do

    end subroutine subblock


end module mod_make_subblock
