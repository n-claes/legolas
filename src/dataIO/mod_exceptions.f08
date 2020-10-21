! =============================================================================
!> Module to explicitly handle exceptions. Depending on the application at hand
!! we override what happens when an exception is raised, which is useful for
!! testing purposes (no error stop if we _expect_ something to fail).
!! Loosely based on an example given in
!! https://github.com/Goddard-Fortran-Ecosystem/pFUnit_demos/blob/main/Basic/src/throw.F90
module mod_exceptions
  implicit none

  private

  abstract interface
    subroutine raise(message)
       character(len=*), intent(in) :: message
    end subroutine raise
  end interface

  !> pointer to method that is used to raise exceptions
  procedure (raise), pointer :: raise_method => null()
  !> logical to check if <tt>raise_method</tt> pointer is assigned
  logical, save :: initialised = .false.

  public :: raise_exception
  public :: set_raise_method

contains


  !> Private subroutine, sets the pointer to the default
  !! method to be used when raising exceptions.
  subroutine initialise_exceptions()
    raise_method => on_exception_raised
    initialised = .true.
  end subroutine initialise_exceptions


  !> Subroutine meant to be publicly called, sets the
  !! routine to be used when raising exceptions.
  !! Calls the initialisation routine if not already done.
  subroutine set_raise_method(method)
    !> subroutine to be used when an exception is raised
    procedure (raise) :: method

    if (.not. initialised) then
       call initialise_exceptions()
    end if
    raise_method => method
  end subroutine set_raise_method


  !> Raises an exception with a given message.
  !! By default, exceptions terminate program execution.
  !! Calls the initialisation routine if not already done.
  subroutine raise_exception(msg)
    !> message to be used when an exception is raised
    character(len=*), intent(in) :: msg

    if (.not. initialised) then
       call initialise_exceptions()
    end if

    call raise_method(msg)
  end subroutine raise_exception


  !> Workflow that is executed by default when
  !! an exception is raised. The argument <tt>message</tt>
  !! is printed to the console and program execution
  !! is terminated.
  subroutine on_exception_raised(msg)
    use mod_global_variables, only: str_len
    use mod_painting, only: paint_string

    !> message to print to the console when exception is raised
    character(len=*), intent(in) :: msg
    character(len=str_len) :: msg_painted

    call paint_string(" ERROR   | " // msg, "red", msg_painted)
    write(*, *) msg_painted
    error stop
  end subroutine on_exception_raised
end module mod_exceptions
