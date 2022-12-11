! =============================================================================
!> Module to build a single quadblock by creating the four comprising subblocks.
!! This is done for every grid interval, for specific integral elements
!! and spline functions.
module mod_make_subblock
  use mod_dims, only: dims_t
  implicit none

  private

  public :: subblock

contains


  !> Builds the quadblock.
  !! The four subblocks are calculated based on the specified
  !! factors, positions in the block, weights and basis functions.
  !! These four subblocks are then added to the quadblock that is
  !! passed along, and the filled quadblock is returned.
  !! @note  For a 2x2 block at index (i, j) in the top-left subblock we have
  !!        <tt> (2i, 2j) </tt> as index of the bottom-right corner of the 2x2 block.
  !!        This is referred to as the "main index" of the 2x2 block.
  !!        The other three corners are then filled by subtracting one from an index.
  !!        We hence use the following main indices for the four different corners
  !!        of the quadblock:
  !!
  !! - top-left:  <tt>(2i, 2j)</tt>
  !! - top-right: <tt>(2i, 2j + dim_subblock)</tt>
  !! - bot-left:  <tt>(2i + dim_subblock, 2j)</tt>
  !! - bot-right: <tt>(2i + dim_subblock, 2j + dim_subblock)</tt> @endnote
  !! @note  The <tt>positions</tt> array contains the actual positions of the block to fill.
  !!        For example, the first index could be <tt>positions(1, :)=(1, 2)</tt>
  !!        which means that the quantity in <tt>factors(1)</tt> goes into
  !!        the 2x2 block at <tt>(1, 2)</tt>. @endnote
  subroutine subblock( &
    quadblock, factors, positions, curr_weight, spline1, spline2, dims &
  )
    use mod_global_variables, only: dp

    !> the quadblock, filled on exit
    complex(dp), intent(inout)  :: quadblock(:, :)
    !> array containing the factors in every 2x2 block
    complex(dp), intent(in)  :: factors(:)
    !> array containing the positions of the blocks to fill
    integer, intent(in)   :: positions(:, :)
    !> the current weight of the Gaussian quadrature
    real(dp), intent(in)  :: curr_weight
    !> the first basis function
    real(dp), intent(in)  :: spline1(4)
    !> the second basis function
    real(dp), intent(in)  :: spline2(4)
    !> dimensions object
    type(dims_t), intent(in) :: dims

    integer :: i
    integer :: curr_position(2)
    integer :: idx(2)
    integer :: dim_subblock
    complex(dp)  :: weighted_factors(size(factors))

    idx(:) = 0
    weighted_factors = factors * curr_weight
    dim_subblock = dims%get_dim_subblock()

    do i = 1, size(weighted_factors)
      curr_position = positions(i, :)

      ! Subblock top-left corner
      idx(:) = 2*curr_position(:)
      quadblock(idx(1)-1, idx(2)-1) = quadblock(idx(1)-1, idx(2)-1) + &
                                      spline1(2) * weighted_factors(i) * spline2(2)
      quadblock(idx(1)-1, idx(2)  ) = quadblock(idx(1)-1, idx(2)  ) + &
                                      spline1(2) * weighted_factors(i) * spline2(4)
      quadblock(idx(1),   idx(2)-1) = quadblock(idx(1)  , idx(2)-1) + &
                                      spline1(4) * weighted_factors(i) * spline2(2)
      quadblock(idx(1),   idx(2)  ) = quadblock(idx(1)  , idx(2)  ) + &
                                      spline1(4) * weighted_factors(i) * spline2(4)

      ! Subblock top-right corner
      idx(:) = [2*curr_position(1), 2*curr_position(2) + dim_subblock]
      quadblock(idx(1)-1, idx(2)-1) = quadblock(idx(1)-1, idx(2)-1) + &
                                      spline1(2) * weighted_factors(i) * spline2(1)
      quadblock(idx(1)-1, idx(2)  ) = quadblock(idx(1)-1, idx(2)  ) + &
                                      spline1(2) * weighted_factors(i) * spline2(3)
      quadblock(idx(1)  , idx(2)-1) = quadblock(idx(1)  , idx(2)-1) + &
                                      spline1(4) * weighted_factors(i) * spline2(1)
      quadblock(idx(1)  , idx(2)  ) = quadblock(idx(1)  , idx(2)  ) + &
                                      spline1(4) * weighted_factors(i) * spline2(3)

      ! Subblock bottom-left corner
      idx(:) = [2*curr_position(1) + dim_subblock, 2*curr_position(2)]
      quadblock(idx(1)-1, idx(2)-1) = quadblock(idx(1)-1, idx(2)-1) + &
                                      spline1(1) * weighted_factors(i) * spline2(2)
      quadblock(idx(1)-1, idx(2)  ) = quadblock(idx(1)-1, idx(2)  ) + &
                                      spline1(1) * weighted_factors(i) * spline2(4)
      quadblock(idx(1)  , idx(2)-1) = quadblock(idx(1)  , idx(2)-1) + &
                                      spline1(3) * weighted_factors(i) * spline2(2)
      quadblock(idx(1)  , idx(2)  ) = quadblock(idx(1)  , idx(2)  ) + &
                                      spline1(3) * weighted_factors(i) * spline2(4)

      ! Subblock bottom-right corner
      idx(:) = 2*curr_position(:) + dim_subblock
      quadblock(idx(1)-1, idx(2)-1) = quadblock(idx(1)-1, idx(2)-1) + &
                                      spline1(1) * weighted_factors(i) * spline2(1)
      quadblock(idx(1)-1, idx(2)  ) = quadblock(idx(1)-1, idx(2)  ) + &
                                      spline1(1) * weighted_factors(i) * spline2(3)
      quadblock(idx(1)  , idx(2)-1) = quadblock(idx(1)  , idx(2)-1) + &
                                      spline1(3) * weighted_factors(i) * spline2(1)
      quadblock(idx(1)  , idx(2)  ) = quadblock(idx(1)  , idx(2)  ) + &
                                      spline1(3) * weighted_factors(i) * spline2(3)
    end do

  end subroutine subblock

end module mod_make_subblock
