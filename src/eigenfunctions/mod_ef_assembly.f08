module mod_ef_assembly
  use mod_global_variables, only: dp, ic, NaN
  use mod_settings, only: settings_t
  use mod_grid, only: grid
  use mod_get_indices, only: get_index
  use mod_logging, only: logger, str
  implicit none

  private

  public :: get_ef_eps
  public :: get_ef_deps
  public :: retransform_eigenfunction
  public :: assemble_eigenfunction

contains

  pure function get_ef_eps(settings, ef_grid) result(ef_eps)
    type(settings_t), intent(in) :: settings
    real(dp), intent(in) :: ef_grid(:)
    real(dp) :: ef_eps(size(ef_grid))
    if (settings%grid%get_geometry() == "Cartesian") then
      ef_eps = 1.0_dp
    else if (settings%grid%get_geometry() == "cylindrical") then
      ef_eps = ef_grid
    else
      ef_eps = NaN
    end if
  end function get_ef_eps


  pure function get_ef_deps(settings) result(ef_deps)
    type(settings_t), intent(in) :: settings
    real(dp) :: ef_deps
    if (settings%grid%get_geometry() == "Cartesian") then
      ef_deps = 0.0_dp
    else if (settings%grid%get_geometry() == "cylindrical") then
      ef_deps = 1.0_dp
    else
      ef_deps = NaN
    end if
  end function get_ef_deps


  function retransform_eigenfunction( &
    ef_name, ef_eps, eigenfunction &
  ) result(ef_transformed)
    character(len=*), intent(in) :: ef_name
    real(dp), intent(in) :: ef_eps(:)
    complex(dp), intent(in) :: eigenfunction(:)
    complex(dp) :: ef_transformed(size(eigenfunction))

    select case(ef_name)
    case("rho", "v3", "T", "a2") ! var -> eps * var
      ef_transformed = eigenfunction / ef_eps
    case("v1") ! v1 -> i * eps * v1
      ef_transformed = eigenfunction / (ef_eps * ic)
    case("v2", "a3") ! var -> var
      ! do nothing
      ef_transformed = eigenfunction
    case("a1") ! a1 -> i*a1
      ef_transformed = eigenfunction / ic
    case default
      call logger%error("wrong eigenfunction name during retransform")
    end select
  end function retransform_eigenfunction


  function assemble_eigenfunction( &
    settings, ef_name, ef_grid, state_vector_index, eigenvector, derivative_order &
  ) result(assembled_ef)
    type(settings_t), intent(in) :: settings
    character(len=*), intent(in) :: ef_name
    real(dp), intent(in) :: ef_grid(:)
    integer, intent(in) :: state_vector_index
    complex(dp), intent(in) :: eigenvector(:)
    complex(dp) :: assembled_ef(settings%grid%get_ef_gridpts())

    integer, intent(in), optional :: derivative_order
    integer :: diff_order, subblock_idx, grid_idx, ef_grid_idx, dim_subblock
    real(dp) :: basis_function(4)

    diff_order = 0
    if (present(derivative_order)) diff_order = derivative_order
    dim_subblock = settings%dims%get_dim_subblock()

    ! map state vector index to actual subblock index, so 1 -> 1, 2 -> 3, 3 -> 5, etc.
    subblock_idx = 2 * state_vector_index - 1

    ! first gridpoint contribution
    basis_function = get_basis_function( &
      ef_name, ef_grid=ef_grid, grid_idx=1, ef_grid_idx=1, diff_order=diff_order &
    )
    assembled_ef(1) = get_combined_value_from_eigenvector( &
      eigenvector, &
      idx=subblock_idx, &
      dim_subblock=dim_subblock, &
      basis_function=basis_function &
    )
    do grid_idx = 1, settings%grid%get_gridpts() - 1
      ! center of interval = 2 * i, end of interval = 2 * i + 1
      do ef_grid_idx = 2 * grid_idx, 2 * grid_idx + 1
        basis_function = get_basis_function( &
          ef_name, &
          ef_grid=ef_grid, &
          grid_idx=grid_idx, &
          ef_grid_idx=ef_grid_idx, &
          diff_order=diff_order &
        )
        assembled_ef(ef_grid_idx) = get_combined_value_from_eigenvector( &
          eigenvector, &
          idx=subblock_idx, &
          dim_subblock=dim_subblock, &
          basis_function=basis_function &
        )
      end do
      subblock_idx = subblock_idx + dim_subblock
    end do
  end function assemble_eigenfunction


  pure complex(dp) function get_combined_value_from_eigenvector( &
    eigenvector, idx, dim_subblock, basis_function &
  )
    complex(dp), intent(in) :: eigenvector(:)
    integer, intent(in) :: idx
    integer, intent(in) :: dim_subblock
    real(dp), intent(in) :: basis_function(:)

    get_combined_value_from_eigenvector = ( &
      eigenvector(idx) * basis_function(2) &
      + eigenvector(idx + 1) * basis_function(4) &
      + eigenvector(idx + dim_subblock) * basis_function(1) &
      + eigenvector(idx + dim_subblock + 1) * basis_function(3) &
    )
  end function get_combined_value_from_eigenvector


  function get_basis_function( &
    ef_name, ef_grid, grid_idx, ef_grid_idx, diff_order &
  ) result(spline)
    use mod_spline_functions

    character(len=*), intent(in) :: ef_name
    real(dp), intent(in) :: ef_grid(:)
    integer, intent(in) :: grid_idx
    integer, intent(in) :: ef_grid_idx
    integer, intent(in) :: diff_order
    real(dp) :: spline(4)
    procedure(), pointer :: basis_function => null()

    select case(ef_name)
    case("v1", "a2", "a3")
      select case(diff_order)
      case(0)
        basis_function => cubic_factors
      case(1)
        basis_function => cubic_factors_deriv
      case(2)
        basis_function => cubic_factors_deriv2
      case default
        call logger%error( &
          "get_basis_function - invalid derivative order given (cubic): " &
          // str(diff_order) &
        )
        return
      end select
    case("rho", "v2", "v3", "T", "a1")
      select case(diff_order)
      case(0)
        basis_function => quadratic_factors
      case(1)
        basis_function => quadratic_factors_deriv
      case default
        call logger%error( &
          "get_basis_function - invalid derivative order given (quadratic): " &
          // str(diff_order) &
        )
        return
      end select
    case default
      call logger%error( &
        "get_basis_function - invalid eigenfunction name given: " // ef_name &
      )
    end select
    call basis_function( &
      ef_grid(ef_grid_idx), grid(grid_idx), grid(grid_idx + 1), spline &
    )
  end function get_basis_function

end module mod_ef_assembly
