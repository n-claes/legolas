module mod_pfunit_exceptions
  use mod_exceptions, only: set_raise_method
  implicit none

  private

  public :: init_pfunit_raise

contains


  subroutine init_pfunit_raise()
    call set_raise_method(on_pfunit_exception_raised)
  end subroutine init_pfunit_raise


  subroutine on_pfunit_exception_raised(message)
    use funit, only: pFUnit_throw => throw

    character(len=*), intent(in) :: message

    call pFUnit_throw(message)
    return
  end subroutine on_pfunit_exception_raised

end module mod_pfunit_exceptions
