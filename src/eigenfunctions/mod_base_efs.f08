module mod_base_efs
  use mod_global_variables, only: dp, str_len_arr
  use mod_settings, only: settings_t
  use mod_get_indices, only: get_index
  use mod_ef_assembly, only: assemble_eigenfunction, get_ef_eps, &
    retransform_eigenfunction
  implicit none

  private

  type, public :: base_ef_t
    character(str_len_arr) :: name
    complex(dp), allocatable :: quantities(:, :)

  contains

    procedure, public :: initialise
    procedure, public :: assemble
    procedure, public :: delete
  end type base_ef_t

contains

  subroutine initialise(this, name, ef_grid_size, nb_efs)
    class(base_ef_t), intent(inout) :: this
    character(str_len_arr), intent(in) :: name
    integer, intent(in) :: ef_grid_size
    integer, intent(in) :: nb_efs

    this%name = name
    allocate(this%quantities(ef_grid_size, nb_efs))
  end subroutine initialise


  subroutine assemble(this, settings, idxs_to_assemble, right_eigenvectors, ef_grid)
    class(base_ef_t), intent(inout) :: this
    type(settings_t), intent(in) :: settings
    integer, intent(in) :: idxs_to_assemble(:)
    complex(dp), intent(in) :: right_eigenvectors(:, :)
    real(dp), intent(in) :: ef_grid(:)

    integer :: i, idx, state_vector_idx
    real(dp), allocatable :: ef_eps(:)
    complex(dp) :: assembled_ef(size(ef_grid))

    ef_eps = get_ef_eps(settings=settings, ef_grid=ef_grid)
    state_vector_idx = get_index(this%name, settings%get_state_vector())
    do i = 1, size(idxs_to_assemble)
      idx = idxs_to_assemble(i)
      assembled_ef = assemble_eigenfunction( &
        settings=settings, &
        ef_name=this%name, &
        ef_grid=ef_grid, &
        state_vector_index=state_vector_idx, &
        eigenvector=right_eigenvectors(:, idx) &
      )
      this%quantities(:, i) = retransform_eigenfunction( &
        ef_name=this%name, &
        ef_eps=ef_eps, &
        eigenfunction=assembled_ef &
      )
    end do
    if (allocated(ef_eps)) deallocate(ef_eps)
  end subroutine assemble


  pure subroutine delete(this)
    class(base_ef_t), intent(inout) :: this
    if(allocated(this%quantities)) deallocate(this%quantities)
  end subroutine delete

end module mod_base_efs
