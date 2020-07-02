! =============================================================================
!> @brief   Module to create the subblocks in a quadblock.
!! @details Builds a single quadblock by creating the four comprising subblocks.
!!          This is done for every grid interval, for specific integral elements
!!          and spline functions.
module mod_make_subblock
  implicit none

  private

  public :: subblock
  public :: reset_factors
  public :: reset_positions


contains


  !> @brief   Builds the quadblock.
  !! @details The four subblocks are calculated based on the specified
  !!          factors, positions in the block, weights and basis functions.
  !!          These four subblocks are then added to the quadblock that is
  !!          passed along, and the filled quadblock is returned.
  !! @note  For a 2x2 block at index (i, j) in the top-left subblock we have
  !!        <tt> (2i, 2j) </tt> as index of the bottom-right corner of the 2x2 block.
  !!        This is referred to as the "main index" of the 2x2 block.
  !!        The other three corners are then filled by subtracting one from an index.
  !!        We hence use the following main indices for the four different corners
  !!        of the quadblock: \n
  !!        - top-left: <tt> (2i, 2j) </tt> \n
  !!        - top-right: <tt> (2i, 2j + dim_subblock) </tt> \n
  !!        - bot-left: <tt> (2i + dim_subblock, 2j) </tt>  \n
  !!        - bot-right: <tt> (2i + dim_subblock, 2j + dim_subblock) </tt>
  !! @note  The \p positions array contains the actual positions of the block to fill.
  !!        For example, the first index could be <tt> positions(1, :) = (1, 2) </tt>
  !!        which means that the quantity in \p factors(1) goes into the 2x2 block at (1, 2).
  !! @param[in, out]  quadblock   the quadblock to fill
  !! @param[in, out]  factors     array containing the factors in each 2x2 block
  !! @param[in] positions         array containing the positions of the block to fill
  !! @param[in] curr_weight       the current weight of the Gaussian quadrature
  !! @param[in] spline1           the first basis function
  !! @param[in] spline2           the second basis function
  subroutine subblock(quadblock, factors, positions, &
                      curr_weight, spline1, spline2)
    use mod_global_variables, only: dp, dim_quadblock, dim_subblock

    complex(dp), intent(inout)  :: quadblock(dim_quadblock, dim_quadblock)
    complex(dp), intent(inout)  :: factors(:)
    integer, intent(in)         :: positions(:, :)
    real(dp), intent(in)        :: curr_weight, spline1(4), spline2(4)

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


  !> @brief   Resets the factors array.
  !! @details If \p factors is allocated, it is deallocated and reallocated
  !!          with size \p size_factors. Afterwards it is initialised to zero.
  !! @param[in, out]  factors   the factors array, reset on exit
  !! @param[in] size_factors    the new size of the factors array
  subroutine reset_factors(factors, size_factors)
    use mod_global_variables, only: dp

    complex(dp), intent(inout), allocatable :: factors(:)
    integer, intent(in)     :: size_factors

    if (allocated(factors)) then
      deallocate(factors)
    end if

    allocate(factors(size_factors))
    factors = (0.0d0, 0.0d0)
  end subroutine reset_factors


  !> @brief   Resets the positions array.
  !! @details If \p positions is allocated, it is deallocated and reallocated
  !!          using \p size_positions. Afterwards it is initialised to zero.
  !!          The new size is always <tt> (size_positions, 2) </tt>.
  !! @param[in, out]  positions   the positions array, reset on exit
  !! @param[in] size_positions    the new size of the positions array
  subroutine reset_positions(positions, size_positions)
    integer, intent(inout), allocatable :: positions(:, :)
    integer, intent(in) :: size_positions

    if (allocated(positions)) then
      deallocate(positions)
    end if

    allocate(positions(size_positions, 2))
    positions = 0
  end subroutine reset_positions

end module mod_make_subblock
