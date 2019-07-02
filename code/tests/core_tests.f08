program core_tests
  use mod_global_variables
  use testmod_assert
  use mod_grid
  use mod_equilibrium
  use mod_equilibrium_derivatives
  use mod_setup_matrix_b
  use mod_setup_matrix_a
  use mod_solvers
  implicit none

  integer             :: test_gridpts

  integer             :: successes
  integer             :: fails
  integer             :: total
  logical             :: bool

  test_gridpts = 4

  call init()

  call test_grid_carth()
  call test_grid_cyl()
  call test_matrix_B()
  call test_matrix_A()
  call test_invert_diagonal_matrix()
  call test_invert_matrix()
  call test_matrix_multiplication()



  call finish()

contains

  !> Initialise test variables.
  subroutine init()
    write(*, *) "=========================="
    write(*, *) "===== RUNNING TESTS ======"
    write(*, *) "=========================="
    !! Set global variables
    call initialise_variables()
    flow = .false.
    radiative_cooling = .false.
    thermal_conduction = .false.
    resistivity = .false.
    external_gravity = .false.

    !! Set grid points
    x_start = 0.0d0
    x_end   = 1.0d0
    call set_gridpts(test_gridpts)

    successes = 0
    fails     = 0
    total     = 0

    write(*, *) "Gridpts         : ", gridpts
    write(*, *) "Gaussian gridpts: ", gauss_gridpts
    write(*, *) "Matrix gridpts  : ", matrix_gridpts
    write(*, *) "x start         : ", int(x_start)
    write(*, *) "x end           : ", int(x_end)
    write(*, *) ""

  end subroutine init


  !> Checks a given test, used for overview and console output.
  subroutine check_test()
    total = total + 1
    if (bool) then
      successes = successes + 1
      write(*, *) "OK"
    else
      fails = fails + 1
      write(*, *) "----> !! FAIL !! <----"
    end if
  end subroutine check_test


  !> Wraps up the test script.
  subroutine finish()
    write(*, *) ""
    write(*, *) "=========================="
    write(*, *) "===== FINISHED TESTS ====="
    write(*, *) "=========================="
    write(*, *) "Total tests     : ", total
    write(*, *) "Tests succeeded : ", successes
    write(*, *) "Tests failed    : ", fails
  end subroutine finish


  !> Tests the Cartesian grid begin and end points
  subroutine test_grid_carth()
    geometry = 'Cartesian'
    call initialise_grid()
    write(*, *) "Testing Cartesian grid start..."
    call assert_real_equal(grid(1), x_start, bool)
    call check_test()
    write(*, *) "Testing Cartesian grid end..."
    call assert_real_equal(grid(gridpts), x_end, bool)
    call check_test()
  end subroutine test_grid_carth


  !> Tests the origin in cylindrical geometry
  subroutine test_grid_cyl()
    geometry = "cylindrical"
    call set_grid_gauss()
    write(*, *) "Testing cylindrical grid start..."
    call assert_real_larger(grid_gauss(1), 0.0d0, bool)
    call check_test()
  end subroutine test_grid_cyl

  !> Test the B-matrix, should be block tri-diagonal
  subroutine test_matrix_B()
    integer     :: i, j, idx_l, idx_r, lb, rb
    real(dp)    :: mat_B(matrix_gridpts, matrix_gridpts)

    write(*, *) "Testing tri-block diagonal matrix B..."

    equilibrium_type = "None"
    call initialise_equilibrium()
    call initialise_equilibrium_derivatives()
    call construct_B(mat_B)

    do i = 1, matrix_gridpts
      do j = 1, matrix_gridpts
        !! This block checks for tri-diagonality:
        !!              lb       rb
        !!  1 <= i <= 16:  1 < j <= 32 not zero
        !! 17 <= i <= 32:  1 < j <= 48 not zero
        !! 33 <= i <= 48: 16 < j <= 64 not zero etc.
        idx_l = (i - 1) / dim_subblock  !! integer division, rounds down
        idx_r = idx_l + 2
        lb = (idx_l - 1) * dim_subblock
        rb = idx_r * dim_subblock
        if (idx_l < 2) then
          lb = 0
        end if
        if (rb > matrix_gridpts) then
          rb = matrix_gridpts
        end if

        if (j <= lb .or. j > rb) then
          call assert_real_equal(mat_B(i, j), 0.0d0, bool)
          if (.not. bool) then
            write(*, *) "    index i, j         : ", i, j
            write(*, *) "    Value of B at index: ", mat_B(i, j)
            write(*, *) "    Value should be    : ", (0.0d0, 0.0d0)
            call check_test()
            call equilibrium_clean()
            call equilibrium_derivatives_clean()
            return
          end if
        end if
      end do
    end do
    call check_test()
    call equilibrium_clean()
    call equilibrium_derivatives_clean()
  end subroutine test_matrix_B


  !> Test the A-matrix, should be block tri-diagonal
  subroutine test_matrix_A()
    integer     :: i, j, idx_l, idx_r, lb, rb
    complex(dp) :: mat_A(matrix_gridpts, matrix_gridpts)

    write(*, *) "Testing tri-block diagonal matrix A..."

    equilibrium_type = "None"
    call initialise_equilibrium()
    call initialise_equilibrium_derivatives()
    call construct_A(mat_A)

    do i = 1, matrix_gridpts
      do j = 1, matrix_gridpts
        !! This block checks for tri-diagonality:
        !!              lb       rb
        !!  1 <= i <= 16:  1 < j <= 32 not zero
        !! 17 <= i <= 32:  1 < j <= 48 not zero
        !! 33 <= i <= 48: 16 < j <= 64 not zero etc.
        idx_l = (i - 1) / dim_subblock  !! integer division, rounds down
        idx_r = idx_l + 2
        lb = (idx_l - 1) * dim_subblock
        rb = idx_r * dim_subblock
        if (idx_l < 2) then
          lb = 1
        end if
        if (rb > matrix_gridpts) then
          rb = matrix_gridpts
        end if

        if (j <= lb .or. j > rb) then
          call assert_complex_equal(mat_A(i, j), (0.0d0, 0.0d0), bool)
          if (.not. bool) then
            write(*, *) "    index i, j         : ", i, j
            write(*, *) "    value of A at index: ", mat_A(i, j)
            write(*, *) "    value should be    : ", (0.0d0, 0.0d0)
            call check_test()
            call equilibrium_clean()
            call equilibrium_derivatives_clean()
            return
          end if
        end if
      end do
    end do
    call check_test()
    call equilibrium_clean()
    call equilibrium_derivatives_clean()
  end subroutine test_matrix_A


  subroutine test_invert_diagonal_matrix()
    real(dp)      :: mat_B(matrix_gridpts, matrix_gridpts)
    real(dp)      :: inv_B(matrix_gridpts, matrix_gridpts), sol_ij
    integer       :: i, j, k

    write(*, *) "Testing inversion of diagonal matrix.."

    mat_B = 0.0d0

    k = 1.0d0
    do i = 1, matrix_gridpts
      do j = 1, matrix_gridpts
        if (i == j) then
          mat_B(i, j) = k
          k = k + 1
        end if
      end do
    end do
    call invert_B(mat_B, inv_B)

    k = 1.0d0
    do i = 1, matrix_gridpts
      do j = 1, matrix_gridpts
        if (i == j) then
          sol_ij = 1.0d0 / k
          k = k + 1
        else
          sol_ij = 0.0d0
        end if
        call assert_real_equal(inv_B(i, j), sol_ij, bool)
        if (.not. bool) then
          write(*, *) "    index i, j          : ", i, j
          write(*, *) "    inverse B calculated: ", inv_B(i, j)
          write(*, *) "    inverse B solution  : ", sol_ij
          call check_test()
          return
        end if
      end do
    end do
    call check_test()
  end subroutine test_invert_diagonal_matrix


  subroutine test_invert_matrix()
    real(dp)    :: mat_B(4, 4)
    real(dp)    :: inv_B(4, 4), inv_B_sol(4, 4)
    integer     :: i, j

    !! Override this just for this purpose, this results in 4x4 matrix
    call set_matrix_gridpts(4)
    write(*, *) "Testing inversion of regular matrix..."
    inv_B = 0.0d0
    mat_B     = reshape((/ 7,  0, -3,  2,   &! column 1
                           2,  3,  4,  2,   &! column 2
                           1, -1, -2, -1,   &
                          -2,  2,  1,  4 /), shape(mat_B))
    inv_B_sol = reshape((/  &
               1.0d0/7.0d0,   0.0d0,        -2.0d0/7.0d0,   -1.0d0/7.0d0,    &
              -4.0d0/7.0d0,   1.0d0,        22.0d0/7.0d0,    4.0d0/7.0d0,    &
              10.0d0/49.0d0, -2.0d0/7.0d0, -76.0d0/49.0d0, -17.0d0/49.0d0,   &
              15.0d0/49.0d0, -3.0d0/7.0d0, -65.0d0/49.0d0, -1.0d0/49.0d0 /), &
                          shape(inv_B_sol))

    call invert_B(mat_B, inv_B)
    do i = 1, 3
      do j = 1, 3
        call assert_real_equal(inv_B(i, j), inv_B_sol(i, j), bool)
        if (.not. bool) then
          write(*, *) "    index i, j          : ", i, j
          write(*, *) "    inverse B calculated: ", inv_B(i, j)
          write(*, *) "    inverse B solution  : ", inv_B_sol(i, j)
          call check_test()
          return
        end if
      end do
    end do
    call check_test()
    !! Reset grid points
    call set_gridpts(test_gridpts)
  end subroutine test_invert_matrix


  subroutine test_matrix_multiplication()
    real(dp)        :: mat_B(4, 4)
    complex(dp)     :: mat_A(4, 4)
    complex(dp)     :: BA(4, 4), BA_sol(4, 4), ir
    integer         :: i, j

    ir = (1.0d0, 0.0d0)

    call set_matrix_gridpts(4)

    write(*, *) "Testing matrix multiplication routine..."
    mat_B = reshape((/ 7.0d0,  0.0d0, -3.0d0,  2.0d0,   &  ! column 1
                       2.0d0,  3.0d0,  4.0d0,  2.0d0,   &  ! column 2
                       1.0d0, -1.0d0, -2.0d0, -1.0d0,   &
                      -2.0d0,  2.0d0,  1.0d0,  4.0d0 /), shape(mat_B))
    mat_A = reshape((/ 2.0d0*ir,  0.0d0*ir,  1.0d0*ic, -3.0d0*ir, &
                       3.0d0*ir,  0.0d0*ir,  2.0d0*ic,  1.0d0*ir, &
                       4.0d0*ic,  3.0d0*ir, -7.0d0*ir,  5.0d0*ic, &
                      -1.0d0*ir, -2.0d0*ir,  3.0d0*ir,  2.0d0*ir /), &
                        shape(mat_A))

    BA_sol = reshape((/ &
      ( 20.0d0, 1.0d0), (-6.0d0,-1.0d0), (-9.0d0,-2.0d0), (-8.0d0,-1.0d0),   &
      ( 19.0d0, 2.0d0), ( 2.0d0,-2.0d0), (-8.0d0,-4.0d0), (10.0d0,-2.0d0),   &
      ( -1.0d0,18.0d0), (16.0d0,10.0d0), (26.0d0,-7.0d0), (13.0d0,28.0d0),   &
      (-12.0d0, 0.0d0), (-5.0d0, 0.0d0), (-9.0d0, 0.0d0), (-1.0d0, 0.0d0)/), &
                      shape(BA_sol))

    call get_B_invA(mat_B, mat_A, BA)

    do i = 1, 4
      do j = 1, 4
        call assert_complex_equal(BA(i, j), BA_sol(i, j), bool)
        if (.not. bool) then
          write(*, *) "    index i, j     : ", i, j
          write(*, *) "    value at index : ", BA(i, j)
          write(*, *) "    value should be: ", BA_sol(i, j)
          call check_test()
          return
        end if
      end do
    end do
    stop
    call check_test()
    call set_gridpts(test_gridpts)
  end subroutine test_matrix_multiplication















end program core_tests
