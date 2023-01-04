module mod_matrix_element_node
  use mod_global_variables, only: dp
  implicit none

  private

  type, public :: matrix_element_node_t
    complex(dp), private :: element
    integer, private :: position(2)
    real(dp), allocatable, private :: spline1(:)
    real(dp), allocatable, private :: spline2(:)
    type(matrix_element_node_t), pointer :: next

  contains

    procedure, public :: get_element
    procedure, public :: get_position
    procedure, public :: get_spline1
    procedure, public :: get_spline2
    procedure, public :: delete
  end type matrix_element_node_t

  public :: new_matrix_element_node

contains

  pure function new_matrix_element_node( &
    element, position, spline1, spline2 &
  ) result(node)
    complex(dp), intent(in) :: element
    integer, intent(in) :: position(2)
    real(dp), intent(in) :: spline1(:)
    real(dp), intent(in) :: spline2(:)
    type(matrix_element_node_t) :: node

    node%element = element
    node%position = position
    allocate(node%spline1, source=spline1)
    allocate(node%spline2, source=spline2)
    node%next => null()
  end function new_matrix_element_node


  pure complex(dp) function get_element(this)
    class(matrix_element_node_t), intent(in) :: this
    get_element = this%element
  end function get_element


  pure function get_position(this) result(position)
    class(matrix_element_node_t), intent(in) :: this
    integer :: position(2)
    position = this%position
  end function get_position


  pure function get_spline1(this) result(spline1)
    class(matrix_element_node_t), intent(in) :: this
    real(dp) :: spline1(size(this%spline1))
    spline1 = this%spline1
  end function get_spline1


  pure function get_spline2(this) result(spline2)
    class(matrix_element_node_t), intent(in) :: this
    real(dp) :: spline2(size(this%spline2))
    spline2 = this%spline2
  end function get_spline2


  pure subroutine delete(this)
    class(matrix_element_node_t), intent(inout) :: this
    if (allocated(this%spline1)) deallocate(this%spline1)
    if (allocated(this%spline2)) deallocate(this%spline2)
    nullify(this%next)
  end subroutine delete

end module mod_matrix_element_node
