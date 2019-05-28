module test_grid
  implicit none

contains

  subroutine init_test_grid(grid, grid_gauss)
    use mod_global_variables
    use mod_setup_grid
    use mod_setup_equilibrium

    double precision, intent(inout)   :: grid(gridpts)
    double precision, intent(out)     :: grid_gauss(4*gridpts)

    call initialise_grid(grid)
    call initialise_equilibrium(grid, grid_gauss)

  end subroutine init_test_grid

end module test_grid
