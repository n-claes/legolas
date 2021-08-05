submodule(mod_eigenfunctions) smod_base_eigenfunctions
  use mod_global_variables, only: ef_gridpts, nb_eqs
  implicit none

contains

  !> Initialises the base eigenfunction array, sets the corresponding names and state
  !! vector indices, allocates the (subset of) eigenfunctions.
  module procedure initialise_base_eigenfunctions
    integer :: i

    ef_names = [ &
      character(len=str_len_arr) :: "rho", "v1", "v2", "v3", "T", "a1", "a2", "a3" &
    ]
    allocate(base_eigenfunctions(size(ef_names)))

    do i = 1, size(base_eigenfunctions)
      base_eigenfunctions(i) % state_vector_index = i
      base_eigenfunctions(i) % name = ef_names(i)
      allocate(base_eigenfunctions(i) % quantities(size(ef_grid), nb_eigenfuncs))
    end do
  end procedure initialise_base_eigenfunctions


  !> Calculates the eigenfunctions corresponding to the requested eigenvalues and
  !! sets them as attributes for the corresponding types.
  module procedure calculate_base_eigenfunctions
    integer :: i, j, eigenvalue_idx
    complex(dp) :: assembled_ef(size(ef_grid))

    do j = 1, size(base_eigenfunctions)
      do i = 1, size(ef_written_idxs)
        eigenvalue_idx = ef_written_idxs(i)
        assembled_ef = get_assembled_eigenfunction( &
          base_eigenfunctions(j), right_eigenvectors(:, eigenvalue_idx) &
        )
        call retransform_eigenfunction(base_eigenfunctions(j) % name, assembled_ef)
        base_eigenfunctions(j) % quantities(:, i) = assembled_ef
      end do
    end do
  end procedure calculate_base_eigenfunctions


  !> Calculates a single eigenfunction based on an eigenvector.
  !! Assembles the eigenfunction corresponding to <tt>ef_idx</tt> and the
  !! given eigenvector. The eigenvector supplied is the one
  !! corresponding to one particular eigenvalue.
  function get_assembled_eigenfunction(base_ef, eigenvector) result(assembled_ef)
    use mod_global_variables, only: dim_subblock, gridpts

    !> the index of the variable in the state vector
    type(ef_type), intent(in)  :: base_ef
    !> the eigenvector for this particular eigenvalue
    complex(dp), intent(in) :: eigenvector(:)
    !> the assembled eigenfunction
    complex(dp) :: assembled_ef(size(ef_grid))
    real(dp)  :: spline(4)
    integer :: subblock_idx, grid_idx, ef_grid_idx

    ! map ef_idx to actual subblock index. So 1 -> 1, 2 -> 3, 3 -> 5 etc.
    subblock_idx = 2 * base_ef%state_vector_index - 1

    ! contribution from first gridpoint
    spline = get_spline(base_ef%name, grid_idx=1, ef_grid_idx=1)
    assembled_ef(1) = get_values_from_eigenvector(eigenvector, subblock_idx, spline)

    do grid_idx = 1, gridpts-1
      ! idx center interval = 2 * i, idx end interval = 2 * i + 1
      do ef_grid_idx = 2 * grid_idx, 2 * grid_idx + 1
        spline = get_spline(base_ef%name, grid_idx, ef_grid_idx)
        assembled_ef(ef_grid_idx) = get_values_from_eigenvector( &
          eigenvector, subblock_idx, spline &
        )
      end do
      ! Increment indices of eigenvector elements
      subblock_idx = subblock_idx + dim_subblock
    end do
  end function get_assembled_eigenfunction


  !> Retrieves the correct values from the eigenvector that correspond to the
  !! requested eigenfunction at the current subblock index mapping.
  function get_values_from_eigenvector(eigenvector, subblock_idx, spline) result(val)
    use mod_global_variables, only: dim_subblock

    !> the eigenvector associated with the requested eigenvalue
    complex(dp), intent(in) :: eigenvector(:)
    !> current subblock index mapping
    integer, intent(in) :: subblock_idx
    !> finite element basis functions corresponding to the current grid position
    real(dp), intent(in)  :: spline(:)
    !> the composed value extracted from the eigenvector
    complex(dp) :: val

    val = ( &
      eigenvector(subblock_idx) * spline(2) &
      + eigenvector(subblock_idx + 1) * spline(4) &
      + eigenvector(subblock_idx + dim_subblock) * spline(1) &
      + eigenvector(subblock_idx + 1 + dim_subblock) * spline(3) &
    )
  end function get_values_from_eigenvector


  !> Returns the finite element basis functions for the given eigenfunction and
  !! position in the grid.
  function get_spline(name, grid_idx, ef_grid_idx) result(spline)
    use mod_spline_functions, only: cubic_factors, quadratic_factors
    use mod_grid, only: grid

    !> name of the current eigenfunction
    character(len=*), intent(in) :: name
    !> current index in the base grid
    integer, intent(in) :: grid_idx
    !> current index in the eigenfunction grid
    integer, intent(in) :: ef_grid_idx
    !> finite element basis functions for this grid position
    real(dp)  :: spline(4)
    procedure(), pointer :: spline_function

    select case(name)
    case("v1", "a2", "a3")
      spline_function => cubic_factors
    case default
      spline_function => quadratic_factors
    end select
    call spline_function( &
      ef_grid(ef_grid_idx), grid(grid_idx), grid(grid_idx + 1), spline &
    )
  end function get_spline


  !> Retransforms the eigenfunctions to their "actual" values.
  !! For example, \(\rho\) is defined as \(\varepsilon\rho\) in the equations,
  !! in this case the eigenfunction is divided by the scale factor.
  subroutine retransform_eigenfunction(name, eigenfunction)
    use mod_global_variables, only: ic
    use mod_logging, only: log_message

    !> name of the current eigenfunction
    character(len=*), intent(in)  :: name
    !> the current eigenfunction, transformed on exit if applicable
    complex(dp), intent(inout)  :: eigenfunction(:)

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
  end subroutine retransform_eigenfunction

end submodule smod_base_eigenfunctions
