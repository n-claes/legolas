module mod_eigenfunctions
  use mod_global_variables, only: dp, str_len_arr, ef_gridpts
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
  !> array with the names of the basis eigenfunctions
  character(str_len_arr), protected, allocatable :: ef_names(:)
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
    module function get_assembled_eigenfunction( &
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
    end function get_assembled_eigenfunction
  end interface

  interface
    module subroutine clean_derived_eigenfunctions(); end subroutine
  end interface

  public :: ef_grid
  public :: ef_names, derived_ef_names
  public :: ef_written_flags, ef_written_idxs
  public :: base_eigenfunctions
  public :: derived_eigenfunctions
  public :: initialise_eigenfunctions
  public :: calculate_eigenfunctions
  public :: eigenfunctions_clean

contains

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


  subroutine calculate_eigenfunctions(right_eigenvectors)
    complex(dp), intent(in) :: right_eigenvectors(:, :)

    call calculate_base_eigenfunctions(right_eigenvectors)
    if (derived_efs_initialised) then
      call calculate_derived_eigenfunctions(right_eigenvectors)
    end if
  end subroutine calculate_eigenfunctions


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


  subroutine eigenfunctions_clean()
    if (efs_initialised) then
      deallocate(ef_grid)
      deallocate(ef_eps)
      deallocate(ef_written_flags)
      deallocate(ef_written_idxs)
      deallocate(ef_names)
      deallocate(base_eigenfunctions)
      if (derived_efs_initialised) then
        call clean_derived_eigenfunctions()
        deallocate(derived_ef_names)
        deallocate(derived_eigenfunctions)
      end if
    end if
    efs_initialised = .false.
    derived_efs_initialised = .false.
  end subroutine eigenfunctions_clean
end module mod_eigenfunctions
