!
! MODULE: mod_make_subblock
!
!> @author
!> Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Module to fill the subblock (quadblock) at every iteration over the
!! Gaussian grid.
!
module mod_make_subblock
  implicit none

  private
  
  public :: subblock
  public :: reset_factors
  public :: reset_positions


contains

  !> Calculates the quadblock at the current iteration over grid_gauss.
  !! @param[in, out] quadblock  The quadblock, filled when returning
  !! @param[in, out] factors    Array containing the integral elements for the
  !!                            current integral parts. Multiplied with current_weight
  !!                            when returning.
  !! @param[in] positions   Array of dimension (2, :) containing the positions
  !!                        of the different integral elements
  !! @param[in] curr_weight Current weight of the Gaussian quadrature
  !! @param[in] spline1   Left basis function of the current integral elements
  !! @param[in] spline2   Right basis function of the current integral elements
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

    !! Remark: positions takes the position of each matrix element in an 8x8
    !! Configuration. This means we have to map the 8x8 positions to the
    !! four 16x16 subblocks in the 32x32 quadblock matrix.
    !! For a main index in the 8x8 matrix given by (i, i), we have:
    !! 1) TOP LEFT SUBBLOCK
    !!    | (2i-1, 2i-1) | (2i-1, 2i) |
    !!    | (2i, 2i-1)   | (2i, 2i)   |
    !! 2) TOP RIGHT SUBBLOCK: use (2i, 2i + dim_subblock) as main indices
    !! 3) BOT LEFT SUBBLOCK : use (2i + dim_subblock, 2i) as main indices
    !! 4) BOT RIGHT SUBBLOCK: use (2i + dim_subblock, 2i + dim_subblock)
    !!                        as main indices
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

  !> Resets the factors array: deallocates the array, reallocates it
  !! with a new size and initialises it to zero.
  !! @param[in, out] factors  The factors array. Out: deallocated and
  !!                          reallocated with size_factors
  !! @param[in] size_factors  The new size of the factors array
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


  !> Resets the positions array: deallocates the array, reallocates it
  !! with a new size and initialises it to zero.
  !! @param[in, out] positions  The positions array, containing the positions
  !!                            of the factors. Out: deallocated and
  !!                            reallocated with size_positions
  !! @param[in] size_positions  The new size of the positions array.
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
