module mod_boundary_manager
  use mod_global_variables, only: dp
  use mod_logging, only: log_message, str
  implicit none

  private

  !> flag to apply essential boundary conditions on T
  logical, save, protected :: apply_T_bounds
  ! flag to apply no-slip boundary conditions on left side (viscosity only)
  logical, save, protected :: apply_noslip_bounds_left
  !> flag to apply no-slip boundary conditions on right side (viscosity only)
  logical, save, protected :: apply_noslip_bounds_right

  interface
    module subroutine apply_essential_boundaries_left(quadblock, matrix)
      complex(dp), intent(inout)  :: quadblock(:, :)
      character, intent(in)  :: matrix
    end subroutine apply_essential_boundaries_left

    module subroutine apply_essential_boundaries_right(quadblock, matrix)
      complex(dp), intent(inout)  :: quadblock(:, :)
      character, intent(in) :: matrix
    end subroutine apply_essential_boundaries_right

    module subroutine apply_natural_boundaries_left(quadblock, matrix)
      complex(dp), intent(inout)  :: quadblock(:, :)
      character, intent(in) :: matrix
    end subroutine apply_natural_boundaries_left

    module subroutine apply_natural_boundaries_right(quadblock, matrix)
      complex(dp), intent(inout)  :: quadblock(:, :)
      character, intent(in) :: matrix
    end subroutine apply_natural_boundaries_right
  end interface

  public :: apply_boundary_conditions

contains

  subroutine apply_boundary_conditions(matrixA, matrixB)
    use mod_global_variables, only: dim_quadblock

    !> the A-matrix with boundary conditions imposed on exit
    complex(dp), intent(inout)  :: matrixA(:, :)
    !> the B-matrix with boundary conditions imposed on exit
    real(dp), intent(inout)     :: matrixB(:, :)
    complex(dp) :: quadblock(dim_quadblock, dim_quadblock)
    integer   :: l_end, r_start

    ! first gridpoint quadblock runs from idx 1 to l_end
    l_end = dim_quadblock
    ! last gridpoint quadblock runs from idx r_start to
    r_start = size(matrixA, dim=1) - dim_quadblock + 1

    call set_boundary_flags()

    ! handle left side boundary conditions B-matrix
    quadblock = matrixB(:l_end, :l_end)
    call apply_natural_boundaries_left(quadblock, matrix="B")
    call apply_essential_boundaries_left(quadblock, matrix="B")
    matrixB(:l_end, :l_end) = real(quadblock)

    ! handle left side boundary conditions A-matrix
    quadblock = matrixA(:l_end, :l_end)
    call apply_natural_boundaries_left(quadblock, matrix="A")
    call apply_essential_boundaries_left(quadblock, matrix="A")
    matrixA(:l_end, :l_end) = quadblock

    ! handle right side boundary conditions B-matrix
    quadblock = matrixB(r_start:, r_start:)
    call apply_natural_boundaries_right(quadblock, matrix="B")
    call apply_essential_boundaries_right(quadblock, matrix="B")
    matrixB(r_start:, r_start:) = real(quadblock)

    ! handle right side boundary conditions A-matrix
    quadblock = matrixA(r_start:, r_start:)
    call apply_natural_boundaries_right(quadblock, matrix="A")
    call apply_essential_boundaries_right(quadblock, matrix="A")
    matrixA(r_start:, r_start:) = quadblock
  end subroutine apply_boundary_conditions


  subroutine set_boundary_flags()
    use mod_equilibrium, only: kappa_field
    use mod_global_variables, only: &
      thermal_conduction, viscosity, coaxial, dp_LIMIT, geometry

    apply_T_bounds = .false.
    apply_noslip_bounds_left = .false.
    apply_noslip_bounds_right = .false.

    ! check if we need regularity conditions on T, this is the case if we have
    ! perpendicular thermal conduction
    if (thermal_conduction .and. any(abs(kappa_field % kappa_perp) > dp_LIMIT)) then
      apply_T_bounds = .true.
    end if

    ! for viscosity, check if we need a no-slip condition.
    if (viscosity) then
      apply_noslip_bounds_right = .true.
      ! does not apply on-axis for cylindrical, unless two coaxial walls are present
      if (coaxial .or. geometry == "Cartesian") then
        apply_noslip_bounds_left = .true.
      end if
    end if
  end subroutine set_boundary_flags

end module mod_boundary_manager
