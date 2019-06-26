!
! MODULE: mod_boundary_conditions
!
!> @author
!> Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Module to calculate the boundary conditions for the eigenvalue problem.
!
module mod_boundary_conditions
  use mod_global_variables
  implicit none

  private

  !> Array containing the 4 quadratic basis functions
  real(dp)    :: h_quadratic(4)
  !> Array containing the derivatives of the 4 quadratic basis functions
  real(dp)    :: dh_quadratic_dr(4)
  !> Array containing the 4 cubic basis functions
  real(dp)    :: h_cubic(4)
  !> Array containing the derivatives of the 4 cubic basis functions
  real(dp)    :: dh_cubic_dr(4)

  public :: boundaries_B_left_edge
  public :: boundaries_B_right_edge
  public :: boundaries_A_left_edge
  public :: boundaries_A_right_edge

contains

  !> Boundary conditions matrix B, left edge (first iteration over quadblock).
  !! B-matrix has no natural boundary conditions, only essential ones.
  !! Dirichlet (fixed) boundary conditions are implemented at the left edge.
  !! @param[in, out] quadblock    The 32x32 quadblock, containing 4 subblocks.
  !!                              Out: boundary conditions applied.
  subroutine boundaries_B_left_edge(quadblock)
    complex(dp), intent(inout) :: quadblock(dim_quadblock, dim_quadblock)

    call fixed_boundaries(quadblock, "l_edge", "B")

  end subroutine boundaries_B_left_edge


  !> Boundary conditions matrix B, right edge (final iteration over quadblock).
  !! B-matrix has no natural boundary conditions, only essential ones.
  !! Dirichlet (fixed) boundary conditions are implemented at the right edge.
  !! @param[in, out] quadblock    The 32x32 quadblock, containing 4 subblocks.
  !!                              Out: boundary conditions applied.
  subroutine boundaries_B_right_edge(quadblock)
    complex(dp), intent(inout) :: quadblock(dim_quadblock, dim_quadblock)

    call fixed_boundaries(quadblock, "r_edge", "B")

  end subroutine boundaries_B_right_edge


  !> Boundary conditions matrix A, left edge. A-matrix has both natural
  !! (from partial integration) and essential (fixed) boundary conditions.
  !! @param[in] eps   The value for epsilon: 1 for Cartesian, r for cylindrical
  !! @param[in, out] quadblock    The 32x32 quadblock, containing 4 subblocks.
  !!                              Out: boundary conditions applied.
  subroutine boundaries_A_left_edge(eps, quadblock)
    real(dp), intent(in)       :: eps
    complex(dp), intent(inout) :: quadblock(dim_quadblock, dim_quadblock)

    call fixed_boundaries(quadblock, "l_edge", "A")

  end subroutine boundaries_A_left_edge


  !> Boundary conditions matrix A, right edge. A-matrix has both natural
  !! (from partial integration) and essential (fixed) boundary conditions.
  !! @param[in] eps   The value for epsilon: 1 for Cartesian, r for cylindrical
  !! @param[in, out] quadblock    The 32x32 quadblock, containing 4 subblocks.
  !!                              Out: boundary conditions applied.
  subroutine boundaries_A_right_edge(eps, quadblock)
    real(dp), intent(in)       :: eps
    complex(dp), intent(inout) :: quadblock(dim_quadblock, dim_quadblock)

    call fixed_boundaries(quadblock, "r_edge", "A")

  end subroutine boundaries_A_right_edge


  !> Boundaries originating from the fixed conditions. When left_edge is
  !! true, set regularity conditions for the left gridpoint. If left_edge is
  !! false, set them for the right gridpoint.
  !! @param[in, out] quadblock    The 32x32 quadblock, containing 4 subblocks.
  !!                              Out: fixed boundary conditions applied.
  !! @param[in] edge    'r_edge' for right boundary, 'l_edge' for left boundary
  !! @param[in] matrix  'A' for A-matrix boundaries, 'B' for B-matrix boundaries
  subroutine fixed_boundaries(quadblock, edge, matrix)
    complex(dp), intent(inout)  :: quadblock(dim_quadblock, dim_quadblock)
    character(6), intent(in)    :: edge
    character, intent(in)       :: matrix
    integer                     :: i
    complex(dp)                 :: unity

    !! \TODO: LEDA sets diagonal elements to ZBIG instead of unity. Are these
    !!        set to zero because of a vary large number??
    if (matrix == "B") then
      unity = (1.0d0, 0.0d0)
    else if (matrix == "A") then
      unity = (0.0d0, 0.0d0)
    else
      write(*, *) "Wrong matrix argument passed to boundaries."
      write(*, *) "Currently set on: ", matrix
      stop
    end if

    !! === LEFT EDGE ===
    !! For the two subbblocks on the left edge of the quadblock, set the
    !! odd columns to zero. For the two subblocks on the top edge of the
    !! quadblock, set the odd rows to zero (such that the right-bottom
    !! subblock is not modified as this one is not at the boundary).
    !! The first diagonal elements of the top-left quadblock are set to unity.
    if (edge == "l_edge") then
      do i = 1, dim_subblock, 2
        ! Every odd column to zero
        quadblock(i, :) = (0.0d0, 0.0d0)
        ! Every odd row to zero
        quadblock(:, i) = (0.0d0, 0.0d0)
        ! Unity on first diagonal elements for B, zero for A
        quadblock(i, i) = unity
      end do
      !! Additional requirements on T1 if thermal conduction is included
      if (thermal_conduction) then
        !! second row + column of T1 set to zero
        quadblock(10, :) = (0.0d0, 0.0d0)
        quadblock(:, 10) = (0.0d0, 0.0d0)
        !! Unity on second diagonal element for B, zero for A
        quadblock(10, 10) = unity
      end if

    !! === RIGHT EDGE ===
    !! For the rigid wall on the right edge, we require v1 = a2 = a3 = 0.
    !! Hence we set the first v1, a1 and a2 columns of the two right subblocks
    !! to zero, and the first v1, a1 and a2 rows of the two bottom subblocks
    !! to zero.
  else if (edge == "r_edge") then
      !! Rows and columns for v1 (index 19 in bottom-right subblock)
      quadblock(19, :) = (0.0d0, 0.0d0)
      quadblock(:, 19) = (0.0d0, 0.0d0)
      !! Insert unity at first diagonal element for B, zero for A
      quadblock(19, 19) = unity

      !! Rows and columns for a2 (index 29 in bottom-right subblock)
      quadblock(29, :) = (0.0d0, 0.0d0)
      quadblock(:, 29) = (0.0d0, 0.0d0)
      !! Insert unity at first diagonal element for B, zero for A
      quadblock(29, 29) = unity

      !! Rows and columns for a3 (index 31 in bottom-right subblock)
      quadblock(31, :) = (0.0d0, 0.0d0)
      quadblock(:, 31) = (0.0d0, 0.0d0)
      !! Insert unity at first diagonal element for B, zero for A
      quadblock(31, 31) = unity

      !! Additional requiremens on T1 if thermal conduction is included
      !! Rows and columns for T1 (2nd) (index 26 in bottom-right subblock)
      quadblock(26, :) = (0.0d0, 0.0d0)
      quadblock(:, 26) = (0.0d0, 0.0d0)
      !! Insert unity at second diagonal element for B, zero for A
      quadblock(26, 26) = unity
    else
      write(*, *) "Wrong edge passed when calling boundaries."
      write(*, *) "Currently set on: ", edge
      stop
    end if

  end subroutine fixed_boundaries

end module mod_boundary_conditions
