program run_tests
  use mod_global_variables
  use mod_assert
  implicit none


  double precision, allocatable :: grid(:), grid_gauss(:)

  !> initialisations
  call init_variables

  allocate(grid(gridpts))
  allocate(grid_gauss(4*gridpts))

  call assert_init
  call tests_init()

  !> test routines
  call test_integral()

  !> get results
  call get_test_results

contains

  subroutine tests_init()
    use test_grid
    call init_test_grid(grid, grid_gauss)
  end subroutine tests_init


  !> simple test routine to check implementation of gaussian weights and nodes
  subroutine test_integral()
    ! take polynomial of degree 7: assume f(x) = 16*x**7
    ! exact integral in [0, 1] is 2
    integer                 :: i
    double precision        :: a, b, sum, xi, solution

    a = 0
    b = 1
    sum = 0
    solution = 2

    do i = 1, n_gauss
      xi = 0.5*(a+b) + 0.5*(b-a)*gaussian_nodes(i)
      sum = sum + 0.5*(b-a)*gaussian_weights(i) * 16*xi**7
    end do

    print*,"routine test_integral()"
    call assert_equal(sum, solution)

  end subroutine test_integral

end program run_tests
