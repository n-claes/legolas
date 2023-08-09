submodule (mod_boundary_manager) smod_natural_boundaries
  use mod_global_variables, only: ic, NaN
  use mod_equilibrium_params, only: k2, k3
  use mod_settings, only: settings_t
  use mod_matrix_elements, only: matrix_elements_t
  use mod_state_vector, only: sv_rho1, sv_v1, sv_v2, sv_v3, sv_T1, sv_a1, sv_a2, sv_a3
  implicit none

  interface
    module subroutine add_natural_regular_terms( &
      x, elements, settings, grid, background &
    )
      real(dp), intent(in) :: x
      type(matrix_elements_t), intent(inout) :: elements
      type(settings_t), intent(in) :: settings
      type(grid_t), intent(in) :: grid
      type(background_t), intent(in) :: background
    end subroutine add_natural_regular_terms

    module subroutine add_natural_flow_terms( &
      x, elements, settings, grid, background &
    )
      real(dp), intent(in) :: x
      type(matrix_elements_t), intent(inout) :: elements
      type(settings_t), intent(in) :: settings
      type(grid_t), intent(in) :: grid
      type(background_t), intent(in) :: background
    end subroutine add_natural_flow_terms

    module subroutine add_natural_resistive_terms( &
      x, elements, settings, grid, background, physics &
    )
      real(dp), intent(in) :: x
      type(matrix_elements_t), intent(inout) :: elements
      type(settings_t), intent(in) :: settings
      type(grid_t), intent(in) :: grid
      type(background_t), intent(in) :: background
      type(physics_t), intent(in) :: physics
    end subroutine add_natural_resistive_terms

    module subroutine add_natural_conduction_terms( &
      x, elements, settings, grid, background, physics &
    )
      real(dp), intent(in) :: x
      type(matrix_elements_t), intent(inout) :: elements
      type(settings_t), intent(in) :: settings
      type(grid_t), intent(in) :: grid
      type(background_t), intent(in) :: background
      type(physics_t), intent(in) :: physics
    end subroutine add_natural_conduction_terms

    module subroutine add_natural_viscosity_terms( &
      x, elements, settings, grid, background &
    )
      real(dp), intent(in) :: x
      type(matrix_elements_t), intent(inout) :: elements
      type(settings_t), intent(in) :: settings
      type(grid_t), intent(in) :: grid
      type(background_t), intent(in) :: background
    end subroutine add_natural_viscosity_terms

    module subroutine add_natural_hall_terms( &
      x, elements, settings, grid, background, physics &
    )
      real(dp), intent(in) :: x
      type(matrix_elements_t), intent(inout) :: elements
      type(settings_t), intent(in) :: settings
      type(grid_t), intent(in) :: grid
      type(background_t), intent(in) :: background
      type(physics_t), intent(in) :: physics
    end subroutine add_natural_hall_terms

    module subroutine add_natural_hall_Bterms( &
      x, elements, settings, grid, background, physics &
    )
      real(dp), intent(in) :: x
      type(matrix_elements_t), intent(inout) :: elements
      type(settings_t), intent(in) :: settings
      type(grid_t), intent(in) :: grid
      type(background_t), intent(in) :: background
      type(physics_t), intent(in) :: physics
    end subroutine add_natural_hall_Bterms
  end interface

contains

  module procedure apply_natural_boundaries_left
    use mod_build_quadblock, only: add_to_quadblock
    use mod_matrix_elements, only: new_matrix_elements

    complex(dp), allocatable :: quadblock(:, :)
    integer :: i, j, dim_quadblock
    real(dp) :: x, weight
    type(matrix_elements_t) :: elements

    dim_quadblock = settings%dims%get_dim_quadblock()
    allocate(quadblock(dim_quadblock, dim_quadblock))
    quadblock = (0.0d0, 0.0d0)

    elements = new_matrix_elements(settings%state_vector)
    x = grid%gaussian_grid(1)
    ! minus one here, since we evaluate boundaries as Bounds[x1] - Bounds[x0]
    weight = -1.0d0

    if (matrix%get_label() == "A") then
      call add_natural_regular_terms(x, elements, settings, grid, background)
      call add_natural_flow_terms(x, elements, settings, grid, background)
      call add_natural_resistive_terms(x, elements, settings, grid, background, physics)
      call add_natural_conduction_terms( &
        x, elements, settings, grid, background, physics &
      )
      call add_natural_viscosity_terms(x, elements, settings, grid, background)
      call add_natural_hall_terms(x, elements, settings, grid, background, physics)
    else if (matrix%get_label() == "B") then
      call add_natural_hall_Bterms(x, elements, settings, grid, background, physics)
    end if
    call add_to_quadblock( &
      quadblock, &
      elements, &
      x, &
      grid%base_grid(1), &
      grid%base_grid(2), &
      weight, &
      settings%dims &
    )

    ! add quadblock elements to left edge
    do j = 1, dim_quadblock
      do i = 1, dim_quadblock
        call matrix%add_element(row=i, column=j, element=quadblock(i, j))
      end do
    end do
    deallocate(quadblock)
    call elements%delete()
  end procedure apply_natural_boundaries_left


  module procedure apply_natural_boundaries_right
    use mod_build_quadblock, only: add_to_quadblock
    use mod_matrix_elements, only: new_matrix_elements

    complex(dp), allocatable :: quadblock(:, :)
    integer :: i, j, ishift, dim_quadblock
    real(dp) :: x, weight
    type(matrix_elements_t) :: elements

    dim_quadblock = settings%dims%get_dim_quadblock()
    allocate(quadblock(dim_quadblock, dim_quadblock))
    quadblock = (0.0d0, 0.0d0)

    elements = new_matrix_elements(settings%state_vector)
    x = grid%gaussian_grid(settings%grid%get_gauss_gridpts())
    weight = 1.0d0

    ! index shift, this is an even number and represents the final index of the
    ! second-to-last quadblock. We add this to the iteration such that it starts
    ! from 1 + ishift, which is an odd number and the starting index of the last
    ! quadblock.
    ishift = matrix%matrix_dim - dim_quadblock

    if (matrix%get_label() == "A") then
      call add_natural_regular_terms(x, elements, settings, grid, background)
      call add_natural_flow_terms(x, elements, settings, grid, background)
      call add_natural_resistive_terms( &
        x, elements, settings, grid, background, physics &
      )
      call add_natural_conduction_terms( &
        x, elements, settings, grid, background, physics &
      )
      call add_natural_viscosity_terms(x, elements, settings, grid, background)
      call add_natural_hall_terms( &
        x, elements, settings, grid, background, physics &
      )
    else if (matrix%get_label() == "B") then
      call add_natural_hall_Bterms( &
        x, elements, settings, grid, background, physics &
      )
    end if
    call add_to_quadblock( &
      quadblock, &
      elements, &
      x, &
      grid%base_grid(settings%grid%get_gridpts() - 1), &
      grid%base_grid(settings%grid%get_gridpts()), &
      weight, &
      settings%dims &
    )

    ! add quadblock elements to right edge
    do j = 1, dim_quadblock
      do i = 1, dim_quadblock
        call matrix%add_element( &
          row=i + ishift, column=j + ishift, element=quadblock(i, j) &
        )
      end do
    end do
    deallocate(quadblock)
    call elements%delete()
  end procedure apply_natural_boundaries_right

end submodule smod_natural_boundaries
