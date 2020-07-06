! =============================================================================
!> Module to build a single quadblock by creating the four comprising subblocks.
!! This is done for every grid interval, for specific integral elements
!! and spline functions.
module mod_make_subblock
  implicit none

  private

  public :: subblock
  public :: reset_factors
  public :: reset_positions


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
  subroutine subblock(quadblock, factors, positions, &
                      curr_weight, spline1, spline2)
    use mod_global_variables, only: dp, dim_quadblock, dim_subblock

    !> the quadblock, filled on exit
    complex(dp), intent(inout)  :: quadblock(dim_quadblock, dim_quadblock)
    !> array containing the factors in every 2x2 block
    complex(dp), intent(inout)  :: factors(:)
    !> array containing the positions of the blocks to fill
    integer, intent(in)         :: positions(:, :)
    !> the current weight of the Gaussian quadrature
    real(dp), intent(in)        :: curr_weight
    !> the first basis function
    real(dp), intent(in)        :: spline1(4)
    !> the second basis function
    real(dp), intent(in)        :: spline2(4)

    integer                     :: i, len_factors
    integer                     :: curr_position(2)
    integer                     :: idx(2)

    idx(:) = 0
    len_factors = size(factors)
    factors = factors * curr_weight

    do i = 1, len_factors
      curr_position = positions(i, :)

      ! Subblock top-left corner
      idx(:) = 2*curr_position(:)
      quadblock(idx(1)-1, idx(2)-1) = quadblock(idx(1)-1, idx(2)-1) + &
                                      spline1(2) * factors(i) * spline2(2)
      quadblock(idx(1)-1, idx(2)  ) = quadblock(idx(1)-1, idx(2)  ) + &
                                      spline1(2) * factors(i) * spline2(4)
      quadblock(idx(1),   idx(2)-1) = quadblock(idx(1)  , idx(2)-1) + &
                                      spline1(4) * factors(i) * spline2(2)
      quadblock(idx(1),   idx(2)  ) = quadblock(idx(1)  , idx(2)  ) + &
                                      spline1(4) * factors(i) * spline2(4)

      ! Subblock top-right corner
      idx(:) = [2*curr_position(1), 2*curr_position(2) + dim_subblock]
      quadblock(idx(1)-1, idx(2)-1) = quadblock(idx(1)-1, idx(2)-1) + &
                                      spline1(2) * factors(i) * spline2(1)
      quadblock(idx(1)-1, idx(2)  ) = quadblock(idx(1)-1, idx(2)  ) + &
                                      spline1(2) * factors(i) * spline2(3)
      quadblock(idx(1)  , idx(2)-1) = quadblock(idx(1)  , idx(2)-1) + &
                                      spline1(4) * factors(i) * spline2(1)
      quadblock(idx(1)  , idx(2)  ) = quadblock(idx(1)  , idx(2)  ) + &
                                      spline1(4) * factors(i) * spline2(3)

      ! Subblock bottom-left corner
      idx(:) = [2*curr_position(1) + dim_subblock, 2*curr_position(2)]
      quadblock(idx(1)-1, idx(2)-1) = quadblock(idx(1)-1, idx(2)-1) + &
                                      spline1(1) * factors(i) * spline2(2)
      quadblock(idx(1)-1, idx(2)  ) = quadblock(idx(1)-1, idx(2)  ) + &
                                      spline1(1) * factors(i) * spline2(4)
      quadblock(idx(1)  , idx(2)-1) = quadblock(idx(1)  , idx(2)-1) + &
                                      spline1(3) * factors(i) * spline2(2)
      quadblock(idx(1)  , idx(2)  ) = quadblock(idx(1)  , idx(2)  ) + &
                                      spline1(3) * factors(i) * spline2(4)

      ! Subblock bottom-right corner
      idx(:) = 2*curr_position(:) + dim_subblock
      quadblock(idx(1)-1, idx(2)-1) = quadblock(idx(1)-1, idx(2)-1) + &
                                      spline1(1) * factors(i) * spline2(1)
      quadblock(idx(1)-1, idx(2)  ) = quadblock(idx(1)-1, idx(2)  ) + &
                                      spline1(1) * factors(i) * spline2(3)
      quadblock(idx(1)  , idx(2)-1) = quadblock(idx(1)  , idx(2)-1) + &
                                      spline1(3) * factors(i) * spline2(1)
      quadblock(idx(1)  , idx(2)  ) = quadblock(idx(1)  , idx(2)  ) + &
                                      spline1(3) * factors(i) * spline2(3)

    end do

  end subroutine subblock


  !> Resets the <tt>factors</tt> array.
  !! If <tt>factors</tt> is allocated, it is deallocated and reallocated
  !! with size <tt>size_factors</tt>. Afterwards it is initialised to zero.
  subroutine reset_factors(factors, size_factors)
    use mod_global_variables, only: dp

    !> the factors array, reset on exit
    complex(dp), intent(inout), allocatable :: factors(:)
    !> the new size of the factors array
    integer, intent(in)     :: size_factors

    if (allocated(factors)) then
      deallocate(factors)
    end if

    allocate(factors(size_factors))
    factors = (0.0d0, 0.0d0)
  end subroutine reset_factors


  !> Resets the <tt>positions</tt> array.
  !! If <tt>positions</tt> is allocated, it is deallocated and reallocated
  !! using <tt>size_positions</tt>. Afterwards it is initialised to zero.
  !! The new size is always <tt> (size_positions, 2) </tt>.
  subroutine reset_positions(positions, size_positions)
    !> the positions array, reset on exit
    integer, intent(inout), allocatable :: positions(:, :)
    !> the new size of the positions array
    integer, intent(in) :: size_positions

    if (allocated(positions)) then
      deallocate(positions)
    end if

    allocate(positions(size_positions, 2))
    positions = 0
  end subroutine reset_positions

end module mod_make_subblock
