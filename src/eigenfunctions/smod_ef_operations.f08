! =============================================================================
!> This submodule contains procedures and functions to be used when we are
!! either assembling eigenfunctions based on the eigenvectors, or retransforming
!! them.
submodule(mod_eigenfunctions) smod_ef_operations
  implicit none

contains

    !> Retransforms the eigenfunctions to their "actual" values.
    !! For example, \(\rho\) is defined as \(\varepsilon\rho\) in the equations,
    !! in this case the eigenfunction is divided by the scale factor.
  module procedure retransform_eigenfunction
    use mod_global_variables, only: ic
    use mod_logging, only: log_message

    select case(name)
    case("rho", "v3", "T", "a2") ! var -> eps * var
      eigenfunction = eigenfunction / ef_eps
    case("v1") ! v1 -> i * eps * v1
      eigenfunction = eigenfunction / (ef_eps * ic)
    case("v2", "a3") ! var -> var
      ! do nothing
    case("a1") ! a1 -> i*a1
      eigenfunction = eigenfunction / ic
    case default
      call log_message("wrong eigenfunction name during retransform", level="error")
    end select
  end procedure retransform_eigenfunction


  !> This function assembles the eigenfunction corresponding to the given
  !! eigenvector. Which eigenfunction to assemble is determined from the
  !! name attribute of the eigenfunction and the derivative order.
  module procedure assemble_eigenfunction
    real(dp)  :: spline(4)
    integer   :: diff_order
    integer :: subblock_idx, grid_idx, ef_grid_idx

    if (present(derivative_order)) then
      diff_order = derivative_order
    else
      diff_order = 0
    end if

    ! map ef_idx to actual subblock index. So 1 -> 1, 2 -> 3, 3 -> 5 etc.
    subblock_idx = 2 * base_ef%state_vector_index - 1

    ! contribution from first gridpoint
    spline = get_spline(base_ef%name, grid_idx=1, ef_grid_idx=1, diff_order=diff_order)
    assembled_ef(1) = get_values_from_eigenvector( &
      eigenvector, subblock_idx, settings%dims, spline &
    )

    do grid_idx = 1, settings%grid%get_gridpts() - 1
      ! idx center interval = 2 * i, idx end interval = 2 * i + 1
      do ef_grid_idx = 2 * grid_idx, 2 * grid_idx + 1
        spline = get_spline(base_ef%name, grid_idx, ef_grid_idx, diff_order=diff_order)
        assembled_ef(ef_grid_idx) = get_values_from_eigenvector( &
          eigenvector, subblock_idx, settings%dims, spline &
        )
      end do
      ! Increment indices of eigenvector elements
      subblock_idx = subblock_idx + settings%dims%get_dim_subblock()
    end do
  end procedure assemble_eigenfunction


  !> Retrieves the correct values from the eigenvector that correspond to the
  !! requested eigenfunction at the current subblock index mapping.
  function get_values_from_eigenvector( &
    eigenvector, subblock_idx, dims, spline &
  ) result(val)
    use mod_dims, only: dims_t

    !> the eigenvector associated with the requested eigenvalue
    complex(dp), intent(in) :: eigenvector(:)
    !> current subblock index mapping
    integer, intent(in) :: subblock_idx
    !> the dimensions object
    type(dims_t), intent(in) :: dims
    !> finite element basis functions corresponding to the current grid position
    real(dp), intent(in)  :: spline(:)
    !> the composed value extracted from the eigenvector
    complex(dp) :: val
    integer :: dim_subblock

    dim_subblock = dims%get_dim_subblock()
    val = ( &
      eigenvector(subblock_idx) * spline(2) &
      + eigenvector(subblock_idx + 1) * spline(4) &
      + eigenvector(subblock_idx + dim_subblock) * spline(1) &
      + eigenvector(subblock_idx + 1 + dim_subblock) * spline(3) &
    )
  end function get_values_from_eigenvector


  !> Returns the finite element basis functions for the given eigenfunction and
  !! position in the grid.
  function get_spline(name, grid_idx, ef_grid_idx, diff_order) result(spline)
    use mod_logging, only: log_message, str
    use mod_spline_functions
    use mod_grid, only: grid

    !> name of the current eigenfunction
    character(len=*), intent(in) :: name
    !> current index in the base grid
    integer, intent(in) :: grid_idx
    !> current index in the eigenfunction grid
    integer, intent(in) :: ef_grid_idx
    !> derivative order of the spline
    integer, intent(in) :: diff_order
    !> finite element basis functions for this grid position
    real(dp)  :: spline(4)
    procedure(), pointer :: spline_function => null()

    select case(name)
    case("v1", "a2", "a3")
      select case(diff_order)
      case(0)
        spline_function => cubic_factors
      case(1)
        spline_function => cubic_factors_deriv
      case(2)
        spline_function => cubic_factors_deriv2
      case default
        call log_message( &
          "get_spline - invalid derivative order given (cubic): " // str(diff_order), &
          level="error" &
        )
        return
      end select

    case("rho", "v2", "v3", "T", "a1")
      select case(diff_order)
      case(0)
        spline_function => quadratic_factors
      case(1)
        spline_function => quadratic_factors_deriv
      case default
        call log_message( &
          "get_spline - invalid order given (quadratic): " // str(diff_order), &
          level="error" &
        )
        return
      end select

    case default
      call log_message( &
        "get_spline - invalid eigenfunction name given: " // name, &
        level="error" &
      )
    end select

    call spline_function( &
      ef_grid(ef_grid_idx), grid(grid_idx), grid(grid_idx + 1), spline &
    )
  end function get_spline

end submodule smod_ef_operations
