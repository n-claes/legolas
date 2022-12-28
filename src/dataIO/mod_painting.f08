! =============================================================================
!> This module handles formatting of terminal-printed strings.
!! Contains subroutines to colourise strings.
module mod_painting
  implicit none

  !> escape character for logging
  character(len=*), parameter :: esc = achar(27)
  !> termination character for logging
  character(len=*), parameter :: term = esc // "[0m"

  !> red ANSI colour sequence
  character(len=*), parameter :: red = esc // "[31m"
  !> green ANSI colour sequence
  character(len=*), parameter :: green = esc // "[32m"
  !> yellow ANSI colour sequence
  character(len=*), parameter :: yellow = esc // "[33m"
  !> blue ANSI colour sequence
  character(len=*), parameter :: blue = esc // "[34m"
  !> cyan ANSI colour sequence
  character(len=*), parameter :: cyan = esc // "[36m"
  !> grey ANSI colour sequence
  character(len=*), parameter :: grey = esc // "[90m"

  private

  public :: paint_string

contains


  ! LCOV_EXCL_START <for convenience, not used when testing>
  !> Subroutine to paint a given string to the desired colour,
  !! returns a new string with ANSI escape sequences prepended
  !! and appended. If the 'colour' argument is not known, simply
  !! returns the string itself.
  pure function paint_string(msg, colour) result(msg_painted)
    !> message to print to the console
    character(len=*), intent(in) :: msg
    !> colour of the message
    character(len=*), intent(in) :: colour
    character(:), allocatable :: fmt
    !> new string with ANSI sequences added
    character(:), allocatable :: msg_painted

    select case(colour)
      case("red")
        fmt = red
      case("green")
        fmt = green
      case("yellow")
        fmt = yellow
      case("blue")
        fmt = blue
      case("cyan")
        fmt = cyan
      case("grey")
        fmt = grey
      case default
        msg_painted = msg
        return
    end select
    msg_painted = trim(adjustl(fmt // msg // term))
  end function paint_string
  ! LCOV_EXCL_STOP
end module mod_painting
