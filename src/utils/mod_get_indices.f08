!=================================================================
!> Module defining convenient index retrieval functions on various arrays.
module mod_get_indices
  implicit none

  private

  !> interface to retrieve the index of an element in an array.
   !! @note: replace get_index by `findloc` once we drop support for gfortran<9 @endnote
  interface get_index
    module procedure find_index_in_character_array
    module procedure find_indices_in_character_array
  end interface get_index


  interface get_subblock_index
    module procedure transform_state_variable_to_subblock_index
  end interface

  public :: get_index
  public :: get_subblock_index

contains

  !> Function to locate the index of a given character in a character array.
  !! Iterates over the elements and returns on the first hit, if no match
  !! was found zero is returned.
  pure function find_index_in_character_array(name, array) result(match_idx)
    !> the name to search for
    character(len=*), intent(in)  :: name
    !> array with the names to search in
    character(len=*), intent(in)  :: array(:)
    !> index of first match
    integer :: match_idx
    integer :: i

    match_idx = 0
    do i = 1, size(array)
      if (array(i) == trim(adjustl(name))) then
        match_idx = i
        exit
      end if
    end do
  end function find_index_in_character_array


  !> Function to locate the indices of an array of characters in another
  !! character array. Returns the indices of the first hit, it no match was
  !! found zero is returned.
  pure function find_indices_in_character_array(names, array) result(match_idxs)
    !> the names to search for
    character(len=*), intent(in)  :: names(:)
    !> array in which to sarch in
    character(len=*), intent(in)  :: array(:)
    !> index of first matches
    integer :: match_idxs(size(names))
    integer :: i

    do i = 1, size(names)
      match_idxs(i) = find_index_in_character_array(names(i), array)
    end do
  end function find_indices_in_character_array


  !> Transforms a given array of state variables to subblock indices.
  !! If `odd` is `.true.` returns the odd indices, otherwise the even indices.
  !! If `edge` equals `"left"`, returns the indices corresponding to the left boundary
  !! block, if `edge` equals `"right"` returns indices for the right boundary block.
  function transform_state_variable_to_subblock_index( &
    variables, state_vector, dim_subblock, odd, edge &
  ) result(idxs)
    use mod_logging, only: logger, str

    !> array of state vector variable names
    character(len=*), intent(in)  :: variables(:)
    !> state vector to consider
    character(len=*), intent(in)  :: state_vector(:)
    !> dimension of the subblock
    integer, intent(in)  :: dim_subblock
    !> uses odd indices if .true. and even if .false.
    logical, intent(in) :: odd
    !> which edge to consider, `"left"` for left boundary and `"right"` for right one
    character(len=*), intent(in)  :: edge
    !> integer array containing the corresponding subblock indices
    integer, allocatable :: idxs(:)

    idxs = get_index(variables, state_vector)
    idxs = 2 * pack(idxs, idxs > 0)

    if (odd) idxs = idxs - 1
    if (edge == "right") idxs = idxs + dim_subblock
    if (size(idxs) == 0) then
      call logger%error( &
        "could not retrieve subblock indices for any variable in " // str(variables) &
        // " for state vector " // str(state_vector) &
      )
      return
    end if
  end function transform_state_variable_to_subblock_index

end module mod_get_indices
