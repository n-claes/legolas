program core_tests
  use mod_global_variables
  implicit none

  call init()

contains

  subroutine init()
    call initialise_variables()
    write(*, *) geometry
    write(*, *) "test ok"

  end subroutine init

end program core_tests
