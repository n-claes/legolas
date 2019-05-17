module mod_global_variables
  implicit none
  public

  character(:), allocatable         ::  geometry
  double precision                  ::  x_start, x_end
  integer                           ::  gridpts, matrix_gridpts


contains

  subroutine init_variables()
    geometry = "cylindrical"
    x_start  = 0.0d0
    x_end    = 1.0d0
    gridpts  = 10
    matrix_gridpts = 16 * gridpts
  end subroutine init_variables

end module mod_global_variables
