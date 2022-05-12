module mod_matrix_node
  use mod_global_variables, only: dp, NaN
  use mod_logging, only: log_message
  implicit none

  private

  !> Base node corresponding to a given column - value pair for a given row
  type, public :: node_t
    !> column index
    integer :: column
    !> value for the matrix element, can be real/complex
    class(*), allocatable :: element
    !> pointer to next node
    class(node_t), pointer :: next

    contains

    procedure :: delete

    procedure :: add_to_node_element
    procedure, private :: add_to_real_node_element
    procedure, private :: add_to_complex_node_element

    procedure :: multiply_node_element_with
    procedure, private :: multiply_real_node_element_with
    procedure, private :: multiply_complex_node_element_with

    procedure, private :: get_real_node_element
    procedure, private :: get_complex_node_element
    generic :: get_node_element => get_real_node_element, get_complex_node_element
  end type node_t

  public :: new_node

contains

  !> Constructor for a new node, sets the column and element attributes.
  !! No nodes are linked yet; the pointer to the next node is initialised to `null()`.
  pure function new_node(column, element) result(node)
    !> column index
    integer, intent(in) :: column
    !> element added to the node
    class(*), intent(in) :: element
    !> new node with given column and element attributes
    type(node_t) :: node

    allocate(node%element, source=element)
    node%column = column
    node%next => null()
  end function new_node


  !> Destructor, deallocates the node attributes.
  pure subroutine delete(this)
    !> type instance
    class(node_t), intent(inout) :: this

    if (allocated(this%element)) deallocate(this%element)
    nullify(this%next)
  end subroutine delete


  !> Generic type-checking function to add a polymorphic element to the existing
  !! element attribute.
  subroutine add_to_node_element(this, element)
    !> type instance
    class(node_t), intent(inout) :: this
    !> element to add to the existing element
    class(*), intent(in) :: element

    select type(node_element => this%element)
      type is (real(dp))
        call this%add_to_real_node_element(element)
      type is (complex(dp))
        call this%add_to_complex_node_element(element)
    end select
  end subroutine add_to_node_element


  !> Sums the current real node element with a new element. The "old" element is
  !! deallocated and reallocated with its new value.
  subroutine add_to_real_node_element(this, element_to_add)
    !> type instance
    class(node_t), intent(inout) :: this
    !> element to add to the existing element
    class(*), intent(in) :: element_to_add
    real(dp) :: existing_element

    call this%get_node_element(existing_element)
    deallocate(this%element)
    select type(element_to_add)
      type is (real(dp))
        allocate(this%element, source=(existing_element + element_to_add))
      type is (complex(dp))
        allocate(this%element, source=(existing_element + element_to_add))
      class default
        allocate(this%element, source=NaN)
        call raise_operation_error(element_type="real", operation="+")
    end select
  end subroutine add_to_real_node_element


  !> Sums the current complex node element with a new element. The "old" element is
  !! deallocated and reallocated with its new value
  subroutine add_to_complex_node_element(this, element_to_add)
    !> type instance
    class(node_t), intent(inout) :: this
    !> element to add to the existing element
    class(*), intent(in) :: element_to_add
    complex(dp) :: existing_element

    call this%get_node_element(existing_element)
    deallocate(this%element)
    select type(element_to_add)
      type is (real(dp))
        allocate(this%element, source=(existing_element + element_to_add))
      type is (complex(dp))
        allocate(this%element, source=(existing_element + element_to_add))
      class default
        allocate(this%element, source=cmplx(NaN, NaN, kind=dp))
        call raise_operation_error(element_type="complex", operation="+")
    end select
  end subroutine add_to_complex_node_element


  !> Generic type-checking routine to multiply a polymorphic node element with
  !! a given number.
  subroutine multiply_node_element_with(this, number)
    !> type instance
    class(node_t), intent(inout) :: this
    !> number to multiply the element with
    class(*), intent(in) :: number

    select type(element => this%element)
      type is (real(dp))
        call this%multiply_real_node_element_with(number)
      type is (complex(dp))
        call this%multiply_complex_node_element_with(number)
    end select
  end subroutine multiply_node_element_with


  !> Multiplies a real node element with a given real/complex number.
  subroutine multiply_real_node_element_with(this, number)
    !> type instance
    class(node_t), intent(inout) :: this
    !> number to multiply the element with
    class(*), intent(in) :: number
    real(dp) :: element

    call this%get_node_element(element)
    deallocate(this%element)
    select type(number)
      type is (real(dp))
        allocate(this%element, source=(number * element))
      type is (complex(dp))
        allocate(this%element, source=(number * element))
      class default
        allocate(this%element, source=NaN)
        call raise_operation_error(element_type="real", operation="*")
    end select
  end subroutine multiply_real_node_element_with


  !> Multiplies a complex node element with a given real/complex number.
  subroutine multiply_complex_node_element_with(this, number)
    !> type instance
    class(node_t), intent(inout) :: this
    !> number to multiply the element with
    class(*), intent(in) :: number
    complex(dp) :: element

    call this%get_node_element(element)
    deallocate(this%element)
    select type(number)
      type is (real(dp))
        allocate(this%element, source=(number * element))
      type is (complex(dp))
        allocate(this%element, source=(number * element))
      class default
        allocate(this%element, source=cmplx(NaN, NaN, kind=dp))
        call raise_operation_error(element_type="complex", operation="*")
    end select
  end subroutine multiply_complex_node_element_with


  !> Getter for nodes with real elements, returns a real `element` attribute.
  !! If the node element is complex it is case to real first, but we throw a
  !! warning since this may indicate unexpected behaviour and information loss.
  !! Throws an error if the element attribute is not of type real/complex.
  subroutine get_real_node_element(this, element)
    !> type instance
    class(node_t), intent(in) :: this
    !> corresponding element
    real(dp), intent(out) :: element

    select type(item => this%element)
      type is (real(dp))
        element = item
      type is (complex(dp))
        element = real(item, kind=dp)
        call raise_type_cast_warning(from_type="complex", to_type="real")
      class default
        element = NaN
        call raise_type_error(element_type="real")
    end select
  end subroutine get_real_node_element


  !> Getter for nodes with complex elements, returns a complex `element` attribute.
  !! If the node element is real instead we cast it to complex, no warning will be
  !! raised since there is no possible loss of information.
  !! Throws an error if the element attribute is not of type real/complex.
  subroutine get_complex_node_element(this, element)
    !> type instance
    class(node_t), intent(in) :: this
    !> corresponding element
    complex(dp), intent(out) :: element

    select type(item => this%element)
      type is (complex(dp))
        element = item
      type is (real(dp))
        element = cmplx(item, kind=dp)
      class default
        element = cmplx(NaN, NaN, kind=dp)
        call raise_type_error(element_type="complex")
    end select
  end subroutine get_complex_node_element


  !> Throws an error message stating that the requested element type and the
  !! corresponding element type do not match.
  subroutine raise_type_error(element_type)
    !> element type to display in the error message
    character(len=*), intent(in) :: element_type

    call log_message("node element is not of type " // element_type, level="error")
  end subroutine raise_type_error


  !> Throws an error message stating that the requested element type and the
  !! corresponding element type do not match.
  subroutine raise_type_cast_warning(from_type, to_type)
    !> original element type to display in the warning message
    character(len=*), intent(in) :: from_type
    !> casted element type
    character(len=*), intent(in) :: to_type

    call log_message( &
      "node element was originally " // from_type &
      // " but has been cast to " // to_type // ". Possible information loss!", &
      level="warning" &
    )
  end subroutine raise_type_cast_warning


  !> Throws an error message stating that element multiplication failed
  subroutine raise_operation_error(element_type, operation)
    !> element type to display in the error message
    character(len=*), intent(in) :: element_type
    !> operation that was performed
    character(len=*), intent(in) :: operation

    call log_message( &
      "unable to do node element (" // element_type // ") " &
      // operation // " given number", &
      level="error" &
    )
  end subroutine raise_operation_error

end module mod_matrix_node
