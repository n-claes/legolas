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

  public

contains

  subroutine boundaries_B_left_edge(quadblock)
    complex(dp), intent(inout) :: quadblock(dim_quadblock, dim_quadblock)

    return

  end subroutine boundaries_B_left_edge

  subroutine boundaries_B_right_edge(quadblock)
    complex(dp), intent(inout) :: quadblock(dim_quadblock, dim_quadblock)

    return

  end subroutine boundaries_B_right_edge

  subroutine boundaries_A_left_edge(quadblock)
    complex(dp), intent(inout) :: quadblock(dim_quadblock, dim_quadblock)

    return

  end subroutine boundaries_A_left_edge

  subroutine boundaries_A_right_edge(quadblock)
    complex(dp), intent(inout) :: quadblock(dim_quadblock, dim_quadblock)

    return

  end subroutine boundaries_A_right_edge

end module mod_boundary_conditions
