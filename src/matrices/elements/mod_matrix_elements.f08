module mod_matrix_elements
  use mod_global_variables, only: dp, NaN
  use mod_matrix_element_node, only: matrix_element_node_t, new_matrix_element_node
  use mod_get_indices, only: get_index
  use mod_logging, only: logger, str
  implicit none

  private

  type, public :: matrix_elements_t
    integer, private :: nb_elements
    integer, private :: spline1_size
    integer, private :: spline2_size
    type(matrix_element_node_t), pointer, private :: head
    type(matrix_element_node_t), pointer, private :: tail
    character(len=:), allocatable, private :: state_vector(:)

  contains

    procedure, public :: add => add_node
    procedure, public :: get_node
    procedure, public :: get_elements
    procedure, public :: get_positions
    procedure, public :: get_nb_elements
    procedure, public :: delete
    procedure, private :: spline_sizes_are_valid
    procedure, private :: increment_nb_elements
  end type matrix_elements_t

  public ::  new_matrix_elements

contains

  pure function new_matrix_elements(state_vector) result(elements)
    character(len=*), intent(in) :: state_vector(:)
    type(matrix_elements_t) :: elements
    elements%nb_elements = 0
    elements%spline1_size = 0
    elements%spline2_size = 0
    elements%head => null()
    elements%tail => null()
    allocate(elements%state_vector, source=state_vector)
  end function new_matrix_elements


  subroutine add_node(this, element, loc1, loc2, spline1, spline2)
    class(matrix_elements_t), intent(inout) :: this
    class(*), intent(in) :: element
    character(len=*), intent(in) :: loc1
    character(len=*), intent(in) :: loc2
    real(dp), intent(in) :: spline1(:)
    real(dp), intent(in) :: spline2(:)
    complex(dp) :: node_element
    integer :: position(2)

    ! position checks
    position(1) = get_index(name=loc1, array=this%state_vector)
    position(2) = get_index(name=loc2, array=this%state_vector)
    if (.not. is_valid_position(position)) return
    ! spline checks
    if (this%spline1_size == 0) this%spline1_size = size(spline1)
    if (this%spline2_size == 0) this%spline2_size = size(spline2)
    if (.not. this%spline_sizes_are_valid(spline1, spline2)) return

    node_element = cast_node_element_to_complex(element)
    ! if head is not associated, then the list is empty and we create the first node
    if (.not. associated(this%head)) then
      allocate( &
        this%head, &
        source=new_matrix_element_node(node_element, position, spline1, spline2) &
      )
      this%tail => this%head
    else
      allocate( &
        this%tail%next, &
        source=new_matrix_element_node(node_element, position, spline1, spline2) &
      )
      this%tail => this%tail%next
    end if
    call this%increment_nb_elements()
  end subroutine add_node


  function get_node(this, inode) result(node)
    class(matrix_elements_t), intent(in) :: this
    integer, intent(in) :: inode
    type(matrix_element_node_t), pointer :: node
    type(matrix_element_node_t), pointer :: current_node
    integer :: i

    node => null()
    if (inode < 1 .or. inode > this%nb_elements) then
      call logger%error("get_node: inode out of range: " // str(inode))
      return
    end if
    current_node => this%head
    do i = 1, inode - 1
      current_node => current_node%next
    end do
    node => current_node
    nullify(current_node)
  end function get_node


  function get_elements(this) result(elements)
    class(matrix_elements_t), intent(in) :: this
    complex(dp) :: elements(this%nb_elements)
    type(matrix_element_node_t), pointer :: current_node
    integer :: i

    current_node => this%head
    do i = 1, this%nb_elements
      elements(i) = current_node%get_element()
      current_node => current_node%next
    end do
    nullify(current_node)
  end function get_elements


  function get_positions(this) result(positions)
    class(matrix_elements_t), intent(in) :: this
    integer :: positions(this%nb_elements, 2)
    type(matrix_element_node_t), pointer :: current_node
    integer :: i

    current_node => this%head
    do i = 1, this%nb_elements
      positions(i, :) = current_node%get_position()
      current_node => current_node%next
    end do
    nullify(current_node)
  end function get_positions


  pure integer function get_nb_elements(this)
    class(matrix_elements_t), intent(in) :: this
    get_nb_elements = this%nb_elements
  end function get_nb_elements


  pure function cast_node_element_to_complex(element) result(node_element)
    class(*), intent(in) :: element
    complex(dp) :: node_element

    select type(element)
      type is (complex(dp))
        node_element = element
      type is (real(dp))
        node_element = cmplx(element, 0.0_dp, kind=dp)
      type is (integer)
        node_element = cmplx(element, 0.0_dp, kind=dp)
      class default
        node_element = cmplx(NaN, NaN, kind=dp)
    end select
  end function cast_node_element_to_complex


  pure logical function is_valid_position(position)
    integer, intent(in) :: position(:)
    is_valid_position = all(position > 0)
  end function is_valid_position


  logical function spline_sizes_are_valid(this, spline1, spline2)
    class(matrix_elements_t), intent(in) :: this
    real(dp), intent(in) :: spline1(:)
    real(dp), intent(in) :: spline2(:)

    spline_sizes_are_valid = .true.
    if (size(spline1) /= this%spline1_size) then
      spline_sizes_are_valid = .false.
      call logger%error( &
        "matrix elements: spline1 size not valid, expected " &
        // str(this%spline1_size) // " but got " // str(size(spline1)) &
      )
    else if (size(spline2) /= this%spline2_size) then
      spline_sizes_are_valid = .false.
      call logger%error( &
        "matrix elements: spline2 size not valid, expected " &
        // str(this%spline2_size) // " but got " // str(size(spline2)) &
      )
    end if
  end function spline_sizes_are_valid


  pure subroutine increment_nb_elements(this)
    class(matrix_elements_t), intent(inout) :: this
    this%nb_elements = this%nb_elements + 1
  end subroutine increment_nb_elements


  pure subroutine delete(this)
    class(matrix_elements_t), intent(inout) :: this
    type(matrix_element_node_t), pointer :: current_node, next_node

    if (.not. associated(this%head)) return

    current_node => this%head
    next_node => current_node%next
    do while (associated(current_node))
      call current_node%delete()
      deallocate(current_node)
      current_node => next_node
      if (associated(current_node)) next_node => current_node%next
    end do
    nullify(current_node)
    nullify(next_node)
    nullify(this%head)
    nullify(this%tail)
    deallocate(this%state_vector)
    this%nb_elements = 0
  end subroutine delete

end module mod_matrix_elements
