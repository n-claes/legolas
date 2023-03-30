module mod_boundary_manager
  use mod_global_variables, only: dp
  use mod_logging, only: logger, str
  use mod_matrix_structure, only: matrix_t
  use mod_matrix_elements, only: matrix_elements_t, new_matrix_elements
  use mod_settings, only: settings_t
  use mod_background, only: background_t
  use mod_physics, only: physics_t
  use mod_grid, only: grid_t
  implicit none

  private

  !> flag to apply essential boundary conditions on T
  logical, save, protected :: apply_T_bounds
  ! flag to apply no-slip boundary conditions on left side (viscosity only)
  logical, save, protected :: apply_noslip_bounds_left
  !> flag to apply no-slip boundary conditions on right side (viscosity only)
  logical, save, protected :: apply_noslip_bounds_right

  interface
    module subroutine apply_essential_boundaries_left(matrix, settings)
      type(matrix_t), intent(inout) :: matrix
      type(settings_t), intent(in) :: settings
    end subroutine apply_essential_boundaries_left

    module subroutine apply_essential_boundaries_right(matrix, settings)
      type(matrix_t), intent(inout) :: matrix
      type(settings_t), intent(in) :: settings
    end subroutine apply_essential_boundaries_right

    module subroutine apply_natural_boundaries_left( &
      matrix, settings, grid, background, physics &
    )
      type(matrix_t), intent(inout) :: matrix
      type(settings_t), intent(in) :: settings
      type(grid_t), intent(in) :: grid
      type(background_t), intent(in) :: background
      type(physics_t), intent(in) :: physics
    end subroutine apply_natural_boundaries_left

    module subroutine apply_natural_boundaries_right( &
      matrix, settings, grid, background, physics &
    )
      type(matrix_t), intent(inout) :: matrix
      type(settings_t), intent(in) :: settings
      type(grid_t), intent(in) :: grid
      type(background_t), intent(in) :: background
      type(physics_t), intent(in) :: physics
    end subroutine apply_natural_boundaries_right
  end interface

  public :: apply_boundary_conditions

contains

  subroutine apply_boundary_conditions( &
    matrix_A, matrix_B, settings, grid, background, physics &
  )
    !> the A-matrix with boundary conditions imposed on exit
    type(matrix_t), intent(inout) :: matrix_A
    !> the B-matrix with boundary conditions imposed on exit
    type(matrix_t), intent(inout) :: matrix_B
    !> the settings object
    type(settings_t), intent(in) :: settings
    type(grid_t), intent(in) :: grid
    type(background_t), intent(in) :: background
    type(physics_t), intent(in) :: physics

    call set_boundary_flags(settings)

    ! handle left side boundary conditions B-matrix
    call apply_natural_boundaries_left(matrix_B, settings, grid, background, physics)
    call apply_essential_boundaries_left(matrix_B, settings)
    ! handle left side boundary conditions A-matrix
    call apply_natural_boundaries_left(matrix_A, settings, grid, background, physics)
    call apply_essential_boundaries_left(matrix_A, settings)
    ! handle right side boundary conditions B-matrix
    call apply_natural_boundaries_right(matrix_B, settings, grid, background, physics)
    call apply_essential_boundaries_right(matrix_B, settings)
    ! handle right side boundary conditions A-matrix
    call apply_natural_boundaries_right(matrix_A, settings, grid, background, physics)
    call apply_essential_boundaries_right(matrix_A, settings)
  end subroutine apply_boundary_conditions


  subroutine set_boundary_flags(settings)
    type(settings_t), intent(in) :: settings


    apply_T_bounds = .false.
    apply_noslip_bounds_left = .false.
    apply_noslip_bounds_right = .false.

    ! check if we need regularity conditions on T, this is the case if we have
    ! perpendicular thermal conduction
    apply_T_bounds = settings%physics%conduction%has_perpendicular_conduction()

    ! for viscosity, check if we need a no-slip condition.
    if (settings%physics%viscosity%is_enabled()) then
      apply_noslip_bounds_right = .true.
      ! does not apply on-axis for cylindrical, unless two coaxial walls are present
      if (settings%grid%coaxial .or. settings%grid%get_geometry() == "Cartesian") then
        apply_noslip_bounds_left = .true.
      end if
    end if
  end subroutine set_boundary_flags

end module mod_boundary_manager
