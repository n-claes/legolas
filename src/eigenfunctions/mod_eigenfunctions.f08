! =============================================================================
!> Main module responsible for eigenfunction management. Contains routines
!! and interfaces to initialise and calculate the eigenfunctions and
!! derived quantities.
module mod_eigenfunctions
  use mod_global_variables, only: dp, str_len_arr, ef_gridpts, state_vector
  use mod_types, only: ef_type
  implicit none

  private

  !> logical to check whether eigenfunctions are initialised
  logical, save, protected :: efs_initialised = .false.
  !> logical to check whether derived eigenfunctions are initialised
  logical, save, protected :: derived_efs_initialised = .false.
  !> grid on which the eigenfunctions are assembled
  real(dp), protected, allocatable :: ef_grid(:)
  !> scale factor dedicated for eigenfunctions (eigenfunction grid != gauss grid)
  real(dp), protected, allocatable :: ef_eps(:)
  !> derivative of scale factor
  real(dp), protected :: ef_deps
  !> array with the names of the derived eigenfunctions
  character(str_len_arr), protected, allocatable :: derived_ef_names(:)
  !> logical array containing flags for written eigenfunctions
  logical, protected, allocatable  :: ef_written_flags(:)
  !> integer array containing indices for written eigenfunctions
  integer, protected, allocatable  :: ef_written_idxs(:)

  !> array with base eigenfunctions
  type(ef_type), allocatable :: base_eigenfunctions(:)
  !> array with derived eigenfunctions
  type(ef_type), allocatable :: derived_eigenfunctions(:)

  interface
    module subroutine initialise_base_eigenfunctions(nb_eigenfuncs)
      integer, intent(in) :: nb_eigenfuncs
    end subroutine initialise_base_eigenfunctions

    module subroutine initialise_derived_eigenfunctions(nb_eigenfuncs)
      integer, intent(in) :: nb_eigenfuncs
    end subroutine initialise_derived_eigenfunctions

    module subroutine calculate_base_eigenfunctions(right_eigenvectors)
      complex(dp), intent(in) :: right_eigenvectors(:, :)
    end subroutine calculate_base_eigenfunctions

    module subroutine calculate_derived_eigenfunctions(right_eigenvectors)
      complex(dp), intent(in) :: right_eigenvectors(:, :)
    end subroutine calculate_derived_eigenfunctions
  end interface

  interface
    module function assemble_eigenfunction( &
      base_ef, eigenvector, derivative_order &
    ) result(assembled_ef)
      !> the base eigenfunction at the current position in the eigenfunction array
      type(ef_type), intent(in) :: base_ef
      !> the eigenvector for the eigenvalue under consideration
      complex(dp), intent(in) :: eigenvector(:)
      !> derivative order of the eigenfunction, defaults to 0
      integer, intent(in), optional :: derivative_order
      !> the assembled eigenfunction (not yet transformed to "actual" values)
      complex(dp) :: assembled_ef(ef_gridpts)
    end function assemble_eigenfunction

    module subroutine retransform_eigenfunction(name, eigenfunction)
      !> name of the current eigenfunction
      character(len=*), intent(in)  :: name
      !> the current eigenfunction, transformed on exit if applicable
      complex(dp), intent(inout)  :: eigenfunction(:)
    end subroutine retransform_eigenfunction
  end interface

  interface
    module subroutine clean_derived_eigenfunctions(); end subroutine
  end interface

  public :: ef_grid
  public :: derived_ef_names
  public :: ef_written_flags, ef_written_idxs
  public :: base_eigenfunctions
  public :: derived_eigenfunctions
  public :: initialise_eigenfunctions
  public :: calculate_eigenfunctions
  public :: retrieve_eigenfunctions
  public :: retrieve_eigenfunction_from_index
  public :: eigenfunctions_clean

contains

  !> Initialises the eigenfunctions based on an array of eigenvalues.
  !! Before initialising all arrays we check which subset of eigenvalues, if any,
  !! needs its eigenfunctions saved.
  subroutine initialise_eigenfunctions(omega)
    use mod_global_variables, only: write_derived_eigenfunctions

    !> the array of calculated eigenvalues
    complex(dp), intent(in) :: omega(:)

    call select_eigenfunctions_to_save(omega)

    call assemble_eigenfunction_grid()
    call initialise_base_eigenfunctions(size(ef_written_idxs))
    if (write_derived_eigenfunctions) then
      call initialise_derived_eigenfunctions(size(ef_written_idxs))
    end if
  end subroutine initialise_eigenfunctions


  !> Calculates both the base eigenfunctions and the derived quantities thereof
  !! based on a 2D array of right eigenvectors.
  subroutine calculate_eigenfunctions(right_eigenvectors)
    complex(dp), intent(in) :: right_eigenvectors(:, :)

    call calculate_base_eigenfunctions(right_eigenvectors)
    if (derived_efs_initialised) then
      call calculate_derived_eigenfunctions(right_eigenvectors)
    end if
  end subroutine calculate_eigenfunctions


  !> Returns the full set of eigenfunctions corresponding to the given eigenfunction
  !! name.
  function retrieve_eigenfunctions(name) result(eigenfunctions)
    use mod_logging, only: log_message
    use mod_get_indices, only: get_index

    !> name of the eigenfunction to retrieve
    character(len=*), intent(in)  :: name
    !> the eigenfunctions corresponding to the given name
    type(ef_type) :: eigenfunctions
    integer :: name_idx

    ! check if we want a regular eigenfunction
    name_idx = get_index(name, state_vector)
    if (name_idx > 0) then
      ! found, retrieve and return
      eigenfunctions = base_eigenfunctions(name_idx)
      return
    end if
    ! not found (= 0), try a derived quantity
    if (derived_efs_initialised) then
      name_idx = get_index(name, derived_ef_names)
      if (name_idx > 0) then
        eigenfunctions = derived_eigenfunctions(name_idx)
        return
      end if
    end if
    ! if still not found then something went wrong
    call log_message( &
      "could not retrieve eigenfunction with name " // name, level="error" &
    )
  end function retrieve_eigenfunctions


  !> Retrieves a single eigenfunction based on its index in the attribute
  !! of the main array. For example, if name equals "rho" and ef_idx equals 2
  !! then this routine returns the quantities attribute evaluated at index 2
  !! for the "rho" eigenfunctions.
  !! @note: this routine is needed apart from <tt>retrieve_eigenfunctions</tt>
  !!        because Fortran does not allow access to a derived type's attribute
  !!        directly through the result of a function call. @endnote
  function retrieve_eigenfunction_from_index(name, ef_index) result(eigenfunction)
    !> name of the eigenfunction to retrieve
    character(len=*), intent(in)  :: name
    !> index of the eigenfunction to retrieve
    integer, intent(in) :: ef_index
    !> the eigenfunction corresponding to the given name and index
    complex(dp) :: eigenfunction(size(ef_grid))
    !> all eigenfunctions corresponding to the given name
    type(ef_type) :: eigenfunctions

    eigenfunctions = retrieve_eigenfunctions(name=name)
    eigenfunction = eigenfunctions%quantities(:, ef_index)
  end function retrieve_eigenfunction_from_index


  !> Allocates and assembles the eigenfunction grid, checks the corresponding
  !! scale factor as well.
  subroutine assemble_eigenfunction_grid()
    use mod_global_variables, only: gridpts, ef_gridpts, geometry
    use mod_grid, only: grid

    integer :: grid_idx

    allocate(ef_grid(ef_gridpts))
    ef_grid = 0.0d0

    ! first gridpoint, left edge
    ef_grid(1) = grid(1)
    ! other gridpoints
    do grid_idx = 1, gridpts - 1
      ! position of center point in grid interval
      ef_grid(2 * grid_idx) = 0.5d0 * (grid(grid_idx) + grid(grid_idx + 1))
      ! position of end point in grid interval
      ef_grid(2 * grid_idx + 1) = grid(grid_idx + 1)
    end do

    allocate(ef_eps(size(ef_grid)))
    if (geometry == "Cartesian") then
      ef_eps = 1.0d0
      ef_deps = 0.0d0
    else
      ef_eps = ef_grid
      ef_deps = 1.0d0
    end if
  end subroutine assemble_eigenfunction_grid


  !> Selects a subset of eigenfunctions to be saved.
  subroutine select_eigenfunctions_to_save(omega)
    use mod_global_variables, only: write_eigenfunction_subset

    !> the array of calculated eigenvalues
    complex(dp), intent(in) :: omega(:)
    integer :: i

    ! check which eigenvalues are inside the given radius
    allocate(ef_written_flags(size(omega)))
    if (write_eigenfunction_subset) then
      ef_written_flags = eigenvalue_is_inside_subset_radius(eigenvalue=omega)
    else
      ef_written_flags = .true.
    end if
    ! extract indices of those eigenvalues
    allocate(ef_written_idxs(count(ef_written_flags)))
    ef_written_idxs = pack([(i, i=1, size(omega))], ef_written_flags)
  end subroutine select_eigenfunctions_to_save


  !> Checks if a specific eigenvalue is within the provided subset radius.
  elemental logical function eigenvalue_is_inside_subset_radius(eigenvalue)
    use mod_global_variables, only: eigenfunction_subset_center, &
      eigenfunction_subset_radius

    !> complex value/array to check
    complex(dp), intent(in) :: eigenvalue
    real(dp)  :: distance_from_subset_center

    distance_from_subset_center = sqrt( &
      (aimag(eigenvalue) - aimag(eigenfunction_subset_center)) ** 2 &
      + (real(eigenvalue) - real(eigenfunction_subset_center)) ** 2 &
    )
    eigenvalue_is_inside_subset_radius = ( &
      distance_from_subset_center <= eigenfunction_subset_radius &
    )
  end function eigenvalue_is_inside_subset_radius


  !> Cleaning routine.
  subroutine eigenfunctions_clean()
    if (efs_initialised) then
      deallocate(base_eigenfunctions)
      deallocate(ef_grid)
      deallocate(ef_eps)
      deallocate(ef_written_flags)
      deallocate(ef_written_idxs)
      if (derived_efs_initialised) then
        call clean_derived_eigenfunctions()
        deallocate(derived_eigenfunctions)
        deallocate(derived_ef_names)
      end if
    end if
    efs_initialised = .false.
    derived_efs_initialised = .false.
  end subroutine eigenfunctions_clean
end module mod_eigenfunctions
