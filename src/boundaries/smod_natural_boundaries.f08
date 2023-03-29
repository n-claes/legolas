submodule (mod_boundary_manager) smod_natural_boundaries
  use mod_global_variables, only: ic, NaN
  use mod_build_quadblock, only: add_to_quadblock
  use mod_equilibrium_params, only: k2, k3
  use mod_matrix_structure, only: matrix_t
  use mod_settings, only: settings_t
  implicit none

  !> quadratic basis functions
  real(dp)  :: h_quad(4)
  !> derivative of quadratic basis functions
  real(dp)  :: dh_quad(4)
  !> cubic basis functions
  real(dp)  :: h_cubic(4)
  !> derivative of cubic basis functions
  real(dp)  :: dh_cubic(4)

  interface
    module subroutine add_natural_regular_terms( &
      x, weight, quadblock, settings, grid, background &
    )
      real(dp), intent(in) :: x
      real(dp), intent(in) :: weight
      complex(dp), intent(inout)  :: quadblock(:, :)
      type(settings_t), intent(in) :: settings
      type(grid_t), intent(in) :: grid
      type(background_t), intent(in) :: background
    end subroutine add_natural_regular_terms

    module subroutine add_natural_flow_terms( &
      x, weight, quadblock, settings, grid, background &
    )
      real(dp), intent(in) :: x
      real(dp), intent(in) :: weight
      complex(dp), intent(inout)  :: quadblock(:, :)
      type(settings_t), intent(in) :: settings
      type(grid_t), intent(in) :: grid
      type(background_t), intent(in) :: background
    end subroutine add_natural_flow_terms

    module subroutine add_natural_resistive_terms( &
      x, weight, quadblock, settings, grid, background, physics &
    )
      real(dp), intent(in) :: x
      real(dp), intent(in) :: weight
      complex(dp), intent(inout)  :: quadblock(:, :)
      type(settings_t), intent(in) :: settings
      type(grid_t), intent(in) :: grid
      type(background_t), intent(in) :: background
      type(physics_t), intent(in) :: physics
    end subroutine add_natural_resistive_terms

    module subroutine add_natural_conduction_terms( &
      x, weight, quadblock, settings, grid, background, physics &
    )
      real(dp), intent(in) :: x
      real(dp), intent(in) :: weight
      complex(dp), intent(inout)  :: quadblock(:, :)
      type(settings_t), intent(in) :: settings
      type(grid_t), intent(in) :: grid
      type(background_t), intent(in) :: background
      type(physics_t), intent(in) :: physics
    end subroutine add_natural_conduction_terms

    module subroutine add_natural_viscosity_terms( &
      x, weight, quadblock, settings, grid, background &
    )
      real(dp), intent(in) :: x
      real(dp), intent(in) :: weight
      complex(dp), intent(inout)  :: quadblock(:, :)
      type(settings_t), intent(in) :: settings
      type(grid_t), intent(in) :: grid
      type(background_t), intent(in) :: background
    end subroutine add_natural_viscosity_terms

    module subroutine add_natural_hall_terms( &
      x, weight, quadblock, settings, grid, background, physics &
    )
      real(dp), intent(in) :: x
      real(dp), intent(in) :: weight
      complex(dp), intent(inout)  :: quadblock(:, :)
      type(settings_t), intent(in) :: settings
      type(grid_t), intent(in) :: grid
      type(background_t), intent(in) :: background
      type(physics_t), intent(in) :: physics
    end subroutine add_natural_hall_terms

    module subroutine add_natural_hall_Bterms( &
      x, weight, quadblock, settings, grid, background, physics &
    )
      real(dp), intent(in) :: x
      real(dp), intent(in) :: weight
      complex(dp), intent(inout)  :: quadblock(:, :)
      type(settings_t), intent(in) :: settings
      type(grid_t), intent(in) :: grid
      type(background_t), intent(in) :: background
      type(physics_t), intent(in) :: physics
    end subroutine add_natural_hall_Bterms
  end interface

contains

  module procedure apply_natural_boundaries_left
    complex(dp), allocatable :: quadblock(:, :)
    integer :: i, j, dim_quadblock
    real(dp) :: x, weight

    dim_quadblock = settings%dims%get_dim_quadblock()
    allocate(quadblock(dim_quadblock, dim_quadblock))
    quadblock = (0.0d0, 0.0d0)
    call set_basis_functions(settings=settings, grid=grid, edge="left")

    x = grid%gaussian_grid(1)
    ! minus one here, since we evaluate boundaries as Bounds[x1] - Bounds[x0]
    weight = -1.0d0

    if (matrix%get_label() == "A") then
      call add_natural_regular_terms(x, weight, quadblock, settings, grid, background)
      call add_natural_flow_terms(x, weight, quadblock, settings, grid, background)
      call add_natural_resistive_terms( &
        x, weight, quadblock, settings, grid, background, physics &
      )
      call add_natural_conduction_terms( &
        x, weight, quadblock, settings, grid, background, physics &
      )
      call add_natural_viscosity_terms(x, weight, quadblock, settings, grid, background)
      call add_natural_hall_terms( &
        x, weight, quadblock, settings, grid, background, physics &
      )
    else if (matrix%get_label() == "B") then
      call add_natural_hall_Bterms( &
        x, weight, quadblock, settings, grid, background, physics &
      )
    end if
    ! add quadblock elements to left edge
    do j = 1, dim_quadblock
      do i = 1, dim_quadblock
        call matrix%add_element(row=i, column=j, element=quadblock(i, j))
      end do
    end do
    deallocate(quadblock)
  end procedure apply_natural_boundaries_left


  module procedure apply_natural_boundaries_right
    complex(dp), allocatable :: quadblock(:, :)
    integer :: i, j, ishift, dim_quadblock
    real(dp) :: x, weight

    dim_quadblock = settings%dims%get_dim_quadblock()
    allocate(quadblock(dim_quadblock, dim_quadblock))
    quadblock = (0.0d0, 0.0d0)
    call set_basis_functions(settings=settings, grid=grid, edge="right")

    x = grid%gaussian_grid(settings%grid%get_gauss_gridpts())
    weight = 1.0d0

    ! index shift, this is an even number and represents the final index of the
    ! second-to-last quadblock. We add this to the iteration such that it starts
    ! from 1 + ishift, which is an odd number and the starting index of the last
    ! quadblock.
    ishift = matrix%matrix_dim - dim_quadblock

    if (matrix%get_label() == "A") then
      call add_natural_regular_terms(x, weight, quadblock, settings, grid, background)
      call add_natural_flow_terms(x, weight, quadblock, settings, grid, background)
      call add_natural_resistive_terms( &
        x, weight, quadblock, settings, grid, background, physics &
      )
      call add_natural_conduction_terms( &
        x, weight, quadblock, settings, grid, background, physics &
      )
      call add_natural_viscosity_terms(x, weight, quadblock, settings, grid, background)
      call add_natural_hall_terms( &
        x, weight, quadblock, settings, grid, background, physics &
      )
    else if (matrix%get_label() == "B") then
      call add_natural_hall_Bterms( &
        x, weight, quadblock, settings, grid, background, physics &
      )
    end if
    ! add quadblock elements to right edge
    do j = 1, dim_quadblock
      do i = 1, dim_quadblock
        call matrix%add_element( &
          row=i + ishift, column=j + ishift, element=quadblock(i, j) &
        )
      end do
    end do
    deallocate(quadblock)
  end procedure apply_natural_boundaries_right


  subroutine set_basis_functions(settings, grid, edge)
    use mod_spline_functions, only: quadratic_factors, quadratic_factors_deriv, &
      cubic_factors, cubic_factors_deriv

    type(settings_t), intent(in) :: settings
    type(grid_t), intent(in) :: grid
    character(len=*), intent(in)  :: edge
    real(dp) :: x_pos, x_left, x_right
    integer :: gridpts

    gridpts = settings%grid%get_gridpts()
    if (edge == "left") then
      x_left = grid%base_grid(1)
      x_right = grid%base_grid(2)
      x_pos = x_left
    else if (edge == "right") then
      x_left = grid%base_grid(gridpts - 1)
      x_right = grid%base_grid(gridpts)
      x_pos = x_right
    else
      call logger%error("natural bounds: invalid edge argument " // edge)
    end if

    ! set the basis functions
    call quadratic_factors(x_pos, x_left, x_right, h_quad)
    call quadratic_factors_deriv(x_pos, x_left, x_right, dh_quad)
    call cubic_factors(x_pos, x_left, x_right, h_cubic)
    call cubic_factors_deriv(x_pos, x_left, x_right, dh_cubic)
  end subroutine

end submodule smod_natural_boundaries
