module testmod_core_tests
  use, intrinsic  ::  ieee_arithmetic
  use mod_global_variables

  use testmod_assert

  implicit none

  integer, parameter  :: test_gridpts = 4

  integer             :: successes
  integer             :: fails
  integer             :: total
  logical             :: test_passed

contains

  !> Initialise test variables.
  subroutine init_core()
    use mod_input

    write(*, *) "==============================="
    write(*, *) "===== RUNNING CORE TESTS ======"
    write(*, *) "==============================="

    ! set up default values
    call read_parfile("")

    successes = 0
    fails     = 0
    total     = 0

  end subroutine init_core



  !> Checks a given test, used for overview and console output.
  subroutine check_test()
    total = total + 1
    if (test_passed) then
      successes = successes + 1
      write(*, *) "OK"
    else
      fails = fails + 1
      write(*, *) "----> !! FAIL !! <----"
    end if
  end subroutine check_test



  !> Wraps up the test script.
  subroutine finish_core()
    use mod_grid, only: grid_clean
    use mod_equilibrium, only: equilibrium_clean
    use mod_equilibrium_derivatives, only: equilibrium_derivatives_clean

    call grid_clean()
    call equilibrium_clean()
    call equilibrium_derivatives_clean()
    write(*, *) ""
    write(*, *) "==============================="
    write(*, *) "===== FINISHED CORE TESTS ====="
    write(*, *) "==============================="
    write(*, *) "Total core tests     : ", total
    write(*, *) "Core tests succeeded : ", successes
    write(*, *) "Core tests failed    : ", fails
    write(*, *) "------------------------------------------------"
    write(*, *) ""
  end subroutine finish_core



  subroutine run_core_tests()
    !! Basic tests, grid etc.
    call print_separation("Testing grids")
    call test_grid_edges()
    call test_grid_structure()
    call test_grid_gauss_structure()

    !! Tests for matrix routines
    call print_separation("Testing matrix routines")
    call test_invert_diagonal_matrix()
    call test_invert_matrix()
    call test_matrix_multiplication()
    call test_matrix_multiplication_blas()
    call test_QR()

    !! Tests for the subblock routines
    call test_subblock()

    !! Tests for the boundary conditions
    call test_boundaries()

    !! Tests the different pre-implemented equilibrium configurations
    call print_separation("Testing different equilibrium configurations")
    call test_equilibria()
  end subroutine run_core_tests



  subroutine print_separation(message)
    character(len=*), intent(in)  :: message

    write(*, *) ""
    write(*, *) "------------------------------------------------"
    write(*, *) message
    write(*, *) "------------------------------------------------"
    write(*, *) ""
  end subroutine print_separation



  !> Tests the Cartesian grid begin and end points
  subroutine test_grid_edges()
    use mod_grid

    geometry = 'Cartesian'
    call initialise_grid()
    write(*, *) "Testing grid endpoints, Cartesian..."
    call assert_equal(grid(1), x_start, test_passed)
    write(*, *) "    grid start..."
    call check_test()
    call assert_equal(grid(gridpts), x_end, test_passed)
    write(*, *) "    grid end..."
    call check_test()
    call grid_clean()

    geometry = 'cylindrical'
    call initialise_grid()
    write(*, *) "Testing grid endpoints, cylindrical..."
    call assert_equal(grid(1), 1.0d-5, test_passed)
    write(*, *) "    grid start..."
    call check_test()
    call assert_equal(grid(gridpts), x_end, test_passed)
    write(*, *) "    grid end..."
    call check_test()
    call grid_clean()
  end subroutine test_grid_edges



  !> Tests structure of the grid
  subroutine test_grid_structure()
    use mod_grid

    real(dp)      :: x1, x2
    integer       :: i

    geometry = "Cartesian"
    call initialise_grid()
    write(*, *) "Testing grid structure, Cartesian..."
    do i = 1, gridpts-1
      x1 = grid(i)
      x2 = grid(i+1)
      call assert_smaller(x1, x2, test_passed)
      if (.not. test_passed) then
        write(*, *) "Grid is not monotone!"
        call check_test()
        call grid_clean()
        return
      end if
    end do
    call grid_clean()
    call check_test()

    geometry = "cylindrical"
    call initialise_grid()
    write(*, *) "Testing grid structure, cylindrical..."
    do i = 1, gridpts-1
      x1 = grid(i)
      x2 = grid(i+1)
      call assert_smaller(x1, x2, test_passed)
      if (.not. test_passed) then
        write(*, *) "Grid is not monotone!"
        call check_test()
        call grid_clean()
        return
      end if
    end do
    call grid_clean()
    call check_test()
  end subroutine test_grid_structure



  subroutine test_grid_gauss_structure()
    use mod_grid

    real(dp)      :: x1, x2
    integer       :: i

    geometry = "Cartesian"
    call initialise_grid()
    write(*, *) "Testing grid_gauss structure, Cartesian..."
    do i = 1, gauss_gridpts-1
      x1 = grid_gauss(i)
      x2 = grid_gauss(i+1)
      call assert_smaller(x1, x2, test_passed)
      if (.not. test_passed) then
        write(*, *) "Grid is not monotone!"
        call check_test()
        call grid_clean()
        return
      end if
    end do
    call grid_clean()
    call check_test()

    geometry = "cylindrical"
    call initialise_grid()
    write(*, *) "Testing grid_gauss structure, cylindrical..."
    do i = 1, gauss_gridpts-1
      x1 = grid_gauss(i)
      x2 = grid_gauss(i+1)
      call assert_smaller(x1, x2, test_passed)
      if (.not. test_passed) then
        write(*, *) "Grid is not monotone!"
        call check_test()
        call grid_clean()
        return
      end if
    end do
    call grid_clean()
    call check_test()
  end subroutine test_grid_gauss_structure



  !> Tests the matrix inversion subroutine, inverts the trivial case
  !! of a diagonal matrix.
  subroutine test_invert_diagonal_matrix()
    use mod_solvers, only: invert_B

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
        call assert_equal(inv_B(i, j), sol_ij, test_passed)
        if (.not. test_passed) then
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



  !> More in-depth test of the matrix inversion routine, inverts a
  !! double precision 4x4 matrix.
  subroutine test_invert_matrix()
    use mod_solvers, only: invert_B

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
        call assert_equal(inv_B(i, j), inv_B_sol(i, j), test_passed)
        if (.not. test_passed) then
          write(*, *) "    index i, j          : ", i, j
          write(*, *) "    inverse B calculated: ", inv_B(i, j)
          write(*, *) "    inverse B solution  : ", inv_B_sol(i, j)
          call check_test()
          call set_gridpts(test_gridpts)
          return
        end if
      end do
    end do
    call check_test()
    !! Reset grid points
    call set_gridpts(test_gridpts)
  end subroutine test_invert_matrix



  !> Tests the matrix multiplication subroutine 'matmul', which is built-in
  !! Fortran. Used to compare results to the ones obtained from BLAS.
  !! Essentially tests BA where B is real and A is complex.
  subroutine test_matrix_multiplication()
    use mod_solvers, only: get_B_invA_matmul

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

    call get_B_invA_matmul(mat_B, mat_A, BA)

    do i = 1, 4
      do j = 1, 4
        call assert_equal(BA(i, j), BA_sol(i, j), test_passed)
        if (.not. test_passed) then
          write(*, *) "    index i, j     : ", i, j
          write(*, *) "    value at index : ", BA(i, j)
          write(*, *) "    value should be: ", BA_sol(i, j)
          call check_test()
          call set_gridpts(test_gridpts)
          return
        end if
      end do
    end do
    call check_test()
    !! Reset grid points
    call set_gridpts(test_gridpts)
  end subroutine test_matrix_multiplication



  !> Tests the BLAS routine to perform matrix multiplication. Essentially
  !! tests BA where B is real and A is complex.
  subroutine test_matrix_multiplication_blas()
    use mod_solvers, only: get_B_invA

    real(dp)        :: mat_B(4, 4)
    complex(dp)     :: mat_A(4, 4)
    complex(dp)     :: BA(4, 4), BA_sol(4, 4), ir
    integer         :: i, j

    ir = (1.0d0, 0.0d0)

    call set_matrix_gridpts(4)

    write(*, *) "Testing matrix multiplication routine, using BLAS..."
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
        call assert_equal(BA(i, j), BA_sol(i, j), test_passed)
        if (.not. test_passed) then
          write(*, *) "    index i, j     : ", i, j
          write(*, *) "    value at index : ", BA(i, j)
          write(*, *) "    value should be: ", BA_sol(i, j)
          call check_test()
          call set_gridpts(test_gridpts)
          return
        end if
      end do
    end do
    call check_test()
    !! Reset grid points
    call set_gridpts(test_gridpts)
  end subroutine test_matrix_multiplication_blas



  !> Tests the QR algorithm from LAPACK. Calculates the eigenvalues of a
  !! 4x4 complex matrix.
  subroutine test_QR()
    use mod_solvers, only: solve_QR

    real(dp)        :: mat_B(4, 4)
    complex(dp)     :: mat_A(4, 4)
    complex(dp)     :: omega(4), omega_sol(4), temp
    complex(dp)     :: vl(4, 4), vr(4, 4)
    integer         :: i, j, minidx
    complex(dp)     :: ir

    write(*, *) "Testing QR algorithm..."

    call set_matrix_gridpts(4)
    ir = (1.0d0, 0.0d0)

    !! Problem is wBX = AX, so use unit matrix for B as inversion
    !! inversion routines are already tested at this point
    do i = 1, 4
      do j = 1, 4
        if (i == j) then
          mat_B(i, j) = 1.0d0
        else
          mat_B(i, j) = 0.0d0
        end if
      end do
    end do

    mat_A = reshape((/ 2.0d0*ir, -1.0d0*ir, 0.0d0*ir,  0.0d0*ir, &
                       9.0d0*ir,  2.0d0*ir, 0.0d0*ir,  0.0d0*ir, &
                       0.0d0*ir,  1.0d0*ir, 3.0d0*ir,  1.0d0*ir, &
                       2.0d0*ir,  0.0d0*ir, 0.0d0*ir, -1.0d0*ir/), &
                          shape(mat_A))

    call solve_QR(mat_A, mat_B, omega, vl, vr)
    omega_sol = (/ (-1.0d0, 0.0d0), (2.0d0, -3.0d0), &
                   ( 2.0d0, 3.0d0), (3.0d0,  0.0d0) /)

    !! sort eigenvalues (selection sort, based on real part)
    do i = 1, size(omega)-1
      minidx = minloc(real(omega(i:)), 1) + i - 1
      if (real(omega(i)) > real(omega(minidx))) then
        temp = omega(i)
        omega(i) = omega(minidx)
        omega(minidx) = temp
      end if
    end do

    do i = 1, 4
      call assert_equal(omega(i), omega_sol(i), test_passed)
      if (.not. test_passed) then
        write(*, *) "    eigenvalue from QR  : ", omega(i)
        write(*, *) "    eigenvalue should be: ", omega_sol(i)
        call check_test()
        return
      end if
    end do
    call check_test()
    call set_gridpts(test_gridpts)
  end subroutine test_QR



  subroutine test_subblock()
    use mod_make_subblock, only: subblock
    complex(dp)   :: quadblock(dim_quadblock, dim_quadblock)
    complex(dp)   :: factors(2), a, b
    real(dp)      :: spline1(4), spline2(4), weight
    integer       :: positions(2, 2), i, idx1(16), idx2(16)

    write(*, *) "Testing subblock creation..."
    a = ( 3.0d0, 1.0d0)
    b = (-1.0d0, 5.0d0)

    quadblock = (0.0d0, 0.0d0)

    !! Choose two factors and positions at random
    factors(1) = a
    positions(1, :) = [4, 3]
    factors(2) = b
    positions(2, :) = [7, 5]

    !! Set weight and splines to unity for testing
    weight = 1.0d0
    spline1 = 1.0d0
    spline2 = 1.0d0
    call subblock(quadblock, factors, positions, weight, spline1, spline2)

    !! First factor
    idx1 = (/ 7, 7, 8, 8,  7,  7,  8,  8, 23, 23, 24, 24, 23, 23, 24, 24 /)
    idx2 = (/ 5, 6, 5, 6, 21, 22, 21, 22,  5,  6,  5,  6, 21, 22, 21, 22 /)
    do i = 1, size(idx1)
      call assert_equal(quadblock(idx1(i), idx2(i)), a, test_passed)
      if (.not. test_passed) then
        write(*, *) "    quadblock position: ", idx1(i), idx2(i)
        write(*, *) "    value at position : ", quadblock(idx1(i), idx2(i))
        write(*, *) "    value should be   : ", a
        call check_test()
        return
      end if
    end do

    !! Second factor
    idx1 = (/ 13, 13, 14, 14, 13, 13, 14, 14, 29, 29, 30, 30, 29, 29, 30, 30 /)
    idx2 = (/  9, 10,  9, 10, 25, 26, 25, 26,  9, 10,  9, 10, 25, 26, 25, 26 /)
    do i = 1, size(idx1)
      call assert_equal(quadblock(idx1(i), idx2(i)), b, test_passed)
      if (.not. test_passed) then
        write(*, *) "    quadblock position: ", idx1(i), idx2(i)
        write(*, *) "    value at position : ", quadblock(idx1(i), idx2(i))
        write(*, *) "    value should be   : ", b
        call check_test()
        return
      end if
    end do
    call check_test()
  end subroutine test_subblock



  subroutine test_boundaries()
    call print_separation("Tests for boundary conditions")
    call test_fixed_boundaries_left()
    call test_fixed_boundaries_right()
  end subroutine test_boundaries



  !> Tests the fixed boundary conditions for the left edge
  subroutine test_fixed_boundaries_left()
    use mod_boundary_conditions, only: fixed_boundaries

    complex(dp)                :: qblock_test(dim_quadblock, dim_quadblock)
    integer                    :: i, j

    write(*, *) "Testing fixed boundary conditions left..."

    qblock_test = (5.0d0, 2.0d0)

    call fixed_boundaries(qblock_test, 'l_edge', 'B')

    do i = 1, dim_subblock, 2
      do j = 1, dim_quadblock
        if (i == j) then
          call assert_equal(qblock_test(i, j), (1.0d0, 0.0d0), test_passed)
        else
          ! Check rows
          call assert_equal(qblock_test(i, j), (0.0d0, 0.0d0), test_passed)
          ! Check columns
          call assert_equal(qblock_test(j, i), (0.0d0, 0.0d0), test_passed)
        end if
        if (.not. test_passed) then
          write(*, *) "    quadblock index  : ", i, j, "/", j, i
          write(*, *) "    value at position: ", qblock_test(i, j), &
                                              "/", qblock_test(j, i)
          if (i == j) then
            write(*, *) "    value should be  : ", (1.0d0, 0.0d0)
          else
            write(*, *) "    value should be  : ", (0.0d0, 0.0d0)
          end if
          call check_test()
          return
        end if
      end do
    end do
    call check_test()
  end subroutine test_fixed_boundaries_left



  !> Tests the fixed boundary conditions for the right edge
  subroutine test_fixed_boundaries_right()
    use mod_boundary_conditions, only: fixed_boundaries

    complex(dp)                :: qblock_test(dim_quadblock, dim_quadblock)
    integer                    :: i, j, curr_idx, idxs(3)

    write(*, *) "Testing fixed boundary conditions right..."
    idxs = (/ 19, 29, 31 /)
    qblock_test = (5.0d0, 2.0d0)
    call fixed_boundaries(qblock_test, 'r_edge', 'B')

    if (boundary_type == 'wall') then

      do i = 1, dim_quadblock
        do j = 1, size(idxs)
          curr_idx = idxs(j)
          if (i == curr_idx) then
            call assert_equal(qblock_test(i, curr_idx), (1.0d0, 0.0d0), &
                                      test_passed)
          else
            call assert_equal(qblock_test(i, curr_idx), (0.0d0, 0.0d0), &
                                      test_passed)
            call assert_equal(qblock_test(curr_idx, i), (0.0d0, 0.0d0), &
                                      test_passed)
          end if
          if (.not. test_passed) then
            write(*, *) "    quadblock index  : ", i, curr_idx, "/", curr_idx, i
            write(*, *) "    value at position: ", qblock_test(i, curr_idx), &
                                                "/", qblock_test(curr_idx, i)
            if (i == curr_idx) then
              write(*, *) "    value should be  : ", (1.0d0, 0.0d0)
            else
              write(*, *) "    value should be  : ", (0.0d0, 0.0d0)
            end if
            call check_test()
            return
          end if
        end do
      end do

    else
      do i = 1, dim_quadblock
        do j = 1, dim_quadblock
          call assert_equal(qblock_test(j, i), (5.0d0, 2.0d0), test_passed)
          if (.not. test_passed) then
            write(*, *) "    quadblock index  : ", j, i
            write(*, *) "    value at position: ", qblock_test(j, i)
            write(*, *) "    value should be  : ", (5.0d0, 2.0d0)
            call check_test()
            return
          end if
        end do
      end do
    end if

    call check_test()
  end subroutine test_fixed_boundaries_right




  !> Tests the different pre-implemented equilibria. Every equilibrium
  !! configuration is checked for :
  !! - Matrices A and B must be block-tri-diagonal
  !! - Matrices A and B can not contain 'inf' or 'NaN' elements
  subroutine test_equilibria()
    character(len=str_len)  :: equil_list(9)
    integer                 :: i

    equil_list = [character(str_len) :: &
                          "Adiabatic homogeneous", &
                          "Resistive homogeneous", &
                          "Gravitational homogeneous", &
                          "Gravito MHD waves", &
                          "Resistive tearing modes", &
                          "Resistive tearing modes with flow", &
                          "Suydam cluster modes", &
                          "Kelvin-Helmholtz", &
                          "Rotating plasma cylinder"]

    do i = 1, size(equil_list)
      equilibrium_type = trim(equil_list(i))
      write(*, *) "EQUILIBRIUM: ", equilibrium_type
      call test_one_equilibrium()
    end do
  end subroutine test_equilibria



  !> Tests one specific equilibrium, this requires that use_precoded is
  !! already set to true, and equilibrium_type is set to one of the
  !! pre-coded cases.
  subroutine test_one_equilibrium()
    use mod_grid, only: initialise_grid, grid_clean
    use mod_equilibrium
    use mod_equilibrium_derivatives
    use mod_matrix_creation, only: create_matrices

    complex(dp)    :: mat_A(matrix_gridpts, matrix_gridpts)
    real(dp)       :: mat_B(matrix_gridpts, matrix_gridpts)

    call initialise_grid()
    call initialise_equilibrium()
    call initialise_equilibrium_derivatives()
    call set_equilibrium()
    call create_matrices(mat_B, mat_A)
    write(*, *) "Testing if matrix B is block-tridiagonal..."
    call test_B_tridiag(mat_B)
    write(*, *) "Testing if matrix B has no 'inf' elements..."
    call test_B_inf(mat_B)
    write(*, *) "Testing if matrix B has no 'NaN' elements..."
    call test_B_nan(mat_B)
    write(*, *) "Testing if matrix B is singular..."
    call test_B_singular(mat_B)
    write(*, *) "Testing if matrix A is block-tridiagonal..."
    call test_A_tridiag(mat_A)
    write(*, *) "Testing if matrix A has no 'inf' elements..."
    call test_A_inf(mat_A)
    write(*, *) "Testing if matrix A has no 'NaN' elements..."
    call test_A_nan(mat_A)
    ! cleanup
    call grid_clean()
    call equilibrium_clean()
    call equilibrium_derivatives_clean()
    write(*, *) "------------------------------------------------"
    write(*, *) ""
  end subroutine test_one_equilibrium



  !> Test the B-matrix, should be block tri-diagonal
  subroutine test_B_tridiag(mat_B)
    real(dp), intent(in)    :: mat_B(matrix_gridpts, matrix_gridpts)
    integer                 :: i, j, idx_l, idx_r, lb, rb

    do i = 1, matrix_gridpts
      do j = 1, matrix_gridpts
        !! This block checks for tri-diagonality:
        !!            lb      rb
        !!  1 <= i <= 16: j > 32 zero
        !! 17 <= i <= 32: j > 48 zero
        !! 33 <= i <= 48: j < 16 and j > 64 zero
        !! 49 <= i <= 64: j < 32 and j > 80 zero etc.
        idx_l = (i - 1) / dim_subblock  !! integer division, rounds down
        idx_r = idx_l + 2
        lb = (idx_l - 1) * dim_subblock + 1
        rb = idx_r * dim_subblock
        if (idx_l < 2) then
          lb = 0
        end if
        if (rb > matrix_gridpts) then
          rb = matrix_gridpts
        end if

        if (j <= lb .or. j > rb) then
          call assert_equal(mat_B(i, j), 0.0d0, test_passed)
          if (.not. test_passed) then
            write(*, *) "    index i, j         : ", i, j
            write(*, *) "    Value of B at index: ", mat_B(i, j)
            write(*, *) "    Value should be    : ", (0.0d0, 0.0d0)
            call check_test()
            return
          end if
        end if
      end do
    end do
    call check_test()
  end subroutine test_B_tridiag


  !> Test the A-matrix, should be block tri-diagonal
  subroutine test_A_tridiag(mat_A)
    complex(dp), intent(in) :: mat_A(matrix_gridpts, matrix_gridpts)
    integer                 :: i, j, idx_l, idx_r, lb, rb

    do i = 1, matrix_gridpts
      do j = 1, matrix_gridpts
        !! This block checks for tri-diagonality:
        !!            lb      rb
        !!  1 <= i <= 16: j > 32 zero
        !! 17 <= i <= 32: j > 48 zero
        !! 33 <= i <= 48: j < 16 and j > 64 zero
        !! 49 <= i <= 64: j < 32 and j > 80 zero etc.
        idx_l = (i - 1) / dim_subblock  !! integer division, rounds down
        idx_r = idx_l + 2
        lb = (idx_l - 1) * dim_subblock + 1
        rb = idx_r * dim_subblock
        if (idx_l < 2) then
          lb = 1
        end if
        if (rb > matrix_gridpts) then
          rb = matrix_gridpts
        end if

        if (j < lb .or. j > rb) then
          call assert_equal(mat_A(i, j), (0.0d0, 0.0d0), test_passed)
          if (.not. test_passed) then
            write(*, *) "    index i, j         : ", i, j
            write(*, *) "    value of A at index: ", mat_A(i, j)
            write(*, *) "    value should be    : ", (0.0d0, 0.0d0)
            call check_test()
            return
          end if
        end if
      end do
    end do
    call check_test()
  end subroutine test_A_tridiag



  !> Tests if an element in the B matrix is infinite.
  subroutine test_B_inf(mat_B)
    real(dp), intent(in)    :: mat_B(matrix_gridpts, matrix_gridpts)
    integer                 :: i, j

    do i = 1, matrix_gridpts
      do j = 1, matrix_gridpts
        call assert_is_finite(mat_B(i, j), test_passed)
        if (.not. test_passed) then
          write(*, *) "    index i, j     : ", i, j
          write(*, *) "    B-value at index is infinite"
          call check_test()
          return
        end if
      end do
    end do
    call check_test()
  end subroutine test_B_inf



  !> Tests if an element in the B matrix is NaN.
  subroutine test_B_nan(mat_B)
    real(dp), intent(in)    :: mat_B(matrix_gridpts, matrix_gridpts)
    integer                 :: i, j

    do i = 1, matrix_gridpts
      do j = 1, matrix_gridpts
        call assert_is_no_nan(mat_B(i, j), test_passed)
        if (.not. test_passed) then
          write(*, *) "    index i, j     : ", i, j
          write(*, *) "    B-value at index is NaN"
          call check_test()
          return
        end if
      end do
    end do
    call check_test()
  end subroutine test_B_nan



  !> Tests if an element in the A matrix is infinite. A is complex, so both
  !! the real and complex parts are tested.
  subroutine test_A_inf(mat_A)
    complex(dp), intent(in)    :: mat_A(matrix_gridpts, matrix_gridpts)
    integer                    :: i, j

    do i = 1, matrix_gridpts
      do j = 1, matrix_gridpts
        call assert_is_finite(mat_A(i, j), test_passed)
        if (.not. test_passed) then
          write(*, *) "    index i, j     : ", i, j
          write(*, *) "    value obtained : ", mat_A(i, j)
          write(*, *) "    A-value at index is infinite"
          call check_test()
          return
        end if
      end do
    end do
    call check_test()
  end subroutine test_A_inf



  !> Tests if an element in the A matrix is NaN. A is complex, so both the
  !! real and complex parts are tested.
  subroutine test_A_nan(mat_A)
    complex(dp), intent(in)    :: mat_A(matrix_gridpts, matrix_gridpts)
    integer                    :: i, j

    do i = 1, matrix_gridpts
      do j = 1, matrix_gridpts
        call assert_is_no_nan(mat_A(i, j), test_passed)
        if (.not. test_passed) then
          write(*, *) "    index i, j     : ", i, j
          write(*, *) "    value obtained : ", mat_A(i, j)
          write(*, *) "    A-value at index is NaN"
          call check_test()
          return
        end if
      end do
    end do
    call check_test()
  end subroutine test_A_nan

  !> Tests if the matrix B is singular (row or column equal to zero).
  !! B is inverted, so this is important
  subroutine test_B_singular(mat_B)
    real(dp), intent(in)    :: mat_B(matrix_gridpts, matrix_gridpts)
    integer                 :: i, j, counter

    logical                 :: equal

    counter = 0
    equal   = .false.

    !! Iterate over columns (check for zero rows)
    do i = 1, matrix_gridpts
      do j = 1, matrix_gridpts
        call assert_equal(mat_B(i, j), 0.0d0, equal)
        if (equal) then
          counter = counter + 1
        end if
        if (counter == matrix_gridpts) then ! this means full row of zeroes
          test_passed = .false. ! hence fail test
          write(*, *) "    Matrix B is singular"
          write(*, *) "    At least one zero row encountered"
          write(*, *) "    Row: ", i
          call check_test()
          return
        end if
      end do
      counter = 0
    end do

    counter = 0
    equal   = .false.

    !! Iterate over rows  (check for zero columns)
    do i = 1, matrix_gridpts
      do j = 1, matrix_gridpts
        call assert_equal(mat_B(j, i), 0.0d0, equal)
        if (equal) then
          counter = counter + 1
        end if
        if (counter == matrix_gridpts) then
          test_passed = .false.  !! fail test
          write(*, *) "    Matrix B is singular"
          write(*, *) "    At least one zero column encountered"
          write(*, *) "    Column: ", i
          call check_test()
          return
        end if
      end do
      counter = 0
    end do

    test_passed = .true.
    call check_test()
  end subroutine test_B_singular

end module testmod_core_tests
