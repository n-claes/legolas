submodule (mod_boundary_manager) smod_natural_boundaries
  use mod_global_variables, only: ic, NaN
  use mod_make_subblock, only: subblock
  use mod_grid, only: grid, eps_grid, d_eps_grid_dr
  use mod_equilibrium, only: rho_field, T_field, B_field
  use mod_equilibrium_params, only: k2, k3
  use mod_matrix_structure, only: matrix_t
  use mod_settings, only: settings_t
  implicit none

  !> current position in the grid
  integer   :: grid_idx
  !> weight for the subblock factors, no integrals here so +1 or -1 depending on side
  real(dp)  :: weight
  !> quadratic basis functions
  real(dp)  :: h_quad(4)
  !> derivative of quadratic basis functions
  real(dp)  :: dh_quad(4)
  !> cubic basis functions
  real(dp)  :: h_cubic(4)
  !> derivative of cubic basis functions
  real(dp)  :: dh_cubic(4)
  !> array of factors for the boundary conditions
  complex(dp), allocatable  :: factors(:)
  !> array of positions, governs place in the left and right quadblocks
  integer, allocatable  :: positions(:, :)

  interface
    module subroutine add_natural_regular_terms(quadblock, settings)
      complex(dp), intent(inout)  :: quadblock(:, :)
      type(settings_t), intent(in) :: settings
    end subroutine add_natural_regular_terms

    module subroutine add_natural_flow_terms(quadblock, settings)
      complex(dp), intent(inout)  :: quadblock(:, :)
      type(settings_t), intent(in) :: settings
    end subroutine add_natural_flow_terms

    module subroutine add_natural_resistive_terms(quadblock, settings)
      complex(dp), intent(inout)  :: quadblock(:, :)
      type(settings_t), intent(in) :: settings
    end subroutine add_natural_resistive_terms

    module subroutine add_natural_conduction_terms(quadblock, settings)
      complex(dp), intent(inout)  :: quadblock(:, :)
      type(settings_t), intent(in) :: settings
    end subroutine add_natural_conduction_terms

    module subroutine add_natural_viscosity_terms(quadblock, settings)
      complex(dp), intent(inout)  :: quadblock(:, :)
      type(settings_t), intent(in) :: settings
    end subroutine add_natural_viscosity_terms

    module subroutine add_natural_hall_terms(quadblock, settings)
      complex(dp), intent(inout)  :: quadblock(:, :)
      type(settings_t), intent(in) :: settings
    end subroutine add_natural_hall_terms

    module subroutine add_natural_hall_Bterms(quadblock, settings)
      complex(dp), intent(inout)  :: quadblock(:, :)
      type(settings_t), intent(in) :: settings
    end subroutine add_natural_hall_Bterms
  end interface

contains

  module procedure apply_natural_boundaries_left
    complex(dp), allocatable :: quadblock(:, :)
    integer :: i, j, dim_quadblock

    dim_quadblock = settings%dims%get_dim_quadblock()
    allocate(quadblock(dim_quadblock, dim_quadblock))
    quadblock = (0.0d0, 0.0d0)
    call set_basis_functions(settings=settings, edge="left")

    if (matrix%get_label() == "A") then
      call add_natural_regular_terms(quadblock, settings)
      call add_natural_flow_terms(quadblock, settings)
      call add_natural_resistive_terms(quadblock, settings)
      call add_natural_conduction_terms(quadblock, settings)
      call add_natural_viscosity_terms(quadblock, settings)
      call add_natural_hall_terms(quadblock, settings)
    else if (matrix%get_label() == "B") then
      call add_natural_hall_Bterms(quadblock, settings)
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

    dim_quadblock = settings%dims%get_dim_quadblock()
    allocate(quadblock(dim_quadblock, dim_quadblock))
    quadblock = (0.0d0, 0.0d0)
    call set_basis_functions(settings=settings, edge="right")

    ! index shift, this is an even number and represents the final index of the
    ! second-to-last quadblock. We add this to the iteration such that it starts
    ! from 1 + ishift, which is an odd number and the starting index of the last
    ! quadblock.
    ishift = matrix%matrix_dim - dim_quadblock

    if (matrix%get_label() == "A") then
      call add_natural_regular_terms(quadblock, settings)
      call add_natural_flow_terms(quadblock, settings)
      call add_natural_resistive_terms(quadblock, settings)
      call add_natural_conduction_terms(quadblock, settings)
      call add_natural_viscosity_terms(quadblock, settings)
      call add_natural_hall_terms(quadblock, settings)
    else if (matrix%get_label() == "B") then
      call add_natural_hall_Bterms(quadblock, settings)
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


  subroutine set_basis_functions(settings, edge)
    use mod_spline_functions, only: quadratic_factors, quadratic_factors_deriv, &
      cubic_factors, cubic_factors_deriv

    type(settings_t), intent(in) :: settings
    character(len=*), intent(in)  :: edge
    real(dp) :: x_pos, x_left, x_right
    integer :: gridpts

    gridpts = settings%grid%get_gridpts()
    if (edge == "left") then
      x_left = grid(1)
      x_right = grid(2)
      x_pos = x_left
      grid_idx = 1
      ! minus one here, since we evaluate boundaries as Bounds[x1] - Bounds[x0]
      weight = -1.0d0
    else if (edge == "right") then
      x_left = grid(gridpts - 1)
      x_right = grid(gridpts)
      x_pos = x_right
      grid_idx = settings%grid%get_gauss_gridpts()
      ! plus one here
      weight = 1.0d0
    else
      call log_message("natural bounds: invalid edge argument " // edge, level="error")
    end if

    ! set the basis functions
    call quadratic_factors(x_pos, x_left, x_right, h_quad)
    call quadratic_factors_deriv(x_pos, x_left, x_right, dh_quad)
    call cubic_factors(x_pos, x_left, x_right, h_cubic)
    call cubic_factors_deriv(x_pos, x_left, x_right, dh_cubic)
  end subroutine

    !> Resets the <tt>factors</tt> and <tt>positions</tt> arrays to a given
  !! new size.
  subroutine reset_factor_positions(new_size)
    !> the new size for the <tt>factors</tt> and <tt>positions</tt> arrays
    integer, intent(in) :: new_size

    if (allocated(factors)) then
      deallocate(factors)
    end if
    allocate(factors(new_size))

    if (allocated(positions)) then
      deallocate(positions)
    end if
    allocate(positions(new_size, 2))
  end subroutine reset_factor_positions

end submodule smod_natural_boundaries
