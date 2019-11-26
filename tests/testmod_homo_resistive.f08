module testmod_homo_resistive
  use mod_global_variables
  use mod_grid
  use mod_equilibrium
  use mod_equilibrium_derivatives
  use mod_setup_matrix_a
  use mod_setup_matrix_b
  use mod_input

  implicit none


  real(dp), allocatable        :: mat_B(:, :)
  complex(dp), allocatable     :: mat_A(:, :)
  complex(dp), allocatable     :: omega(:)
  complex(dp), allocatable     :: vl(:, :)
  complex(dp), allocatable     :: vr(:, :)

contains

  subroutine init_homo_resistive_test()
    call read_parfile("tests/testfile_homo_resistive.par")

    call initialise_grid()
    call initialise_equilibrium()
    call initialise_equilibrium_derivatives()
    call set_equilibrium()

    allocate(mat_B(matrix_gridpts, matrix_gridpts))
    allocate(mat_A(matrix_gridpts, matrix_gridpts))
    allocate(omega(matrix_gridpts))
    allocate(vl(matrix_gridpts, matrix_gridpts))
    allocate(vr(matrix_gridpts, matrix_gridpts))
  end subroutine init_homo_resistive_test

  subroutine run_homo_resistive_test()
    use mod_solvers
    use mod_physical_constants
    use mod_io

    complex(dp)   :: mat_A_theory(8, 8)
    real(dp)      :: mat_B_theory(8, 8)

    complex(dp)   :: omega_theory(8)
    complex(dp)   :: vl_theory(8, 8), vr_theory(8, 8)
    real(dp)      :: k1
    real(dp)      :: T0, B03, rho0, eta
    integer       :: n, it_end, i

    real(dp)      :: w_real, w_imag
    logical       :: append

    !! Solve using the code itself
    call construct_B(mat_B)
    call construct_A(mat_A)
    call solve_QR(mat_A, mat_B, omega, vl, vr)

    call save_eigenvalues(omega, "tests/homo_resistive_test_CODE", &
                          append=.false., stream=.true.)

    !! Rely on analytical solution and quantization of wave number,
    !! kx will be quantized according to kx = n*pi / L
    mat_A_theory = 0.0d0
    mat_B_theory = 0.0d0

    ! Set variables (use 1 for index but does not matter as it is homogeneous)
    rho0 = rho0_eq(1)
    T0   = T0_eq(1)
    B03  = B03_eq(1)
    eta  = fixed_eta_value

    ! construct matrix B
    mat_B_theory(1, 1) = 1.0d0
    mat_B_theory(2, 2) = rho0
    mat_B_theory(3, 3) = rho0
    mat_B_theory(4, 4) = rho0
    mat_B_theory(5, 5) = rho0
    mat_B_theory(6, 6) = 1.0d0
    mat_B_theory(7, 7) = 1.0d0
    mat_B_theory(8, 8) = 1.0d0

    ! Code has #eigenvalues = #matrix_gridpts, and every iteration here
    ! yields 8 eigenvalues hence iterate up to it_end
    it_end = matrix_gridpts / 8

    !! To use QR in this case, set gridpoints to 8
    call set_matrix_gridpts(8)

    ! do not append first time
    append = .false.

    do n = 1, it_end
      k1 = n * dpi / (x_end - x_start)

      ! construct matrix A
      mat_A_theory(1, 2) = k1 * rho0
      mat_A_theory(1, 4) = k3 * rho0
      mat_A_theory(2, 1) = T0 * k1
      mat_A_theory(2, 5) = rho0*k1
      mat_A_theory(2, 7) = ic * B03 * (2*k1**2 + k3**2)
      mat_A_theory(3, 6) = -ic * B03 * k3**2
      mat_A_theory(3, 8) = -ic * B03 * k1 * k3
      mat_A_theory(4, 1) = k3 * T0
      mat_A_theory(4, 5) = k3 * rho0
      mat_A_theory(4, 7) = ic * B03 * k1 * k3
      mat_A_theory(5, 2) = gamma_1 * rho0 * T0 * k1
      mat_A_theory(5, 4) = gamma_1 * rho0 * T0 * k3
      mat_A_theory(5, 5) = (0.0d0, 0.0d0) ! resistivity derivative
      mat_A_theory(5, 7) = 2*gamma_1*B03*k1*(k1**2 + k3**2)*eta
      mat_A_theory(6, 3) = ic * B03
      mat_A_theory(6, 6) = -ic * eta * k3**2
      mat_A_theory(6, 8) = ic * eta * k1 * k3
      mat_A_theory(7, 2) = -ic * B03
      mat_A_theory(7, 5) = (0.0d0, 0.0d0) ! resistivity derivative
      mat_A_theory(7, 7) = -ic * eta * (k1**2 + k3**2)
      mat_A_theory(8, 6) = ic * eta * k1 * k3
      mat_A_theory(8, 8) = -ic * eta * k1**2

      ! make it complex
      mat_A_theory = mat_A_theory * (1.0d0, 0.0d0)
      call solve_QR(mat_A_theory, mat_B_theory, omega_theory, vl_theory, &
                    vr_theory)

      !after writing once, start appending to file
      if (n == 2) then
        append = .true.
      end if

      call save_eigenvalues(omega_theory, "tests/homo_resistive_test_THEORY", &
                            append=append, stream=.true.)

    end do

  end subroutine run_homo_resistive_test

  subroutine finish_homo_resistive_test()
    deallocate(mat_B)
    deallocate(mat_A)
    deallocate(omega)
    deallocate(vl)
    deallocate(vr)
    write(*, *) ">> Resistive homogeneous test completed"
    call grid_clean()
    call equilibrium_clean()
    call equilibrium_derivatives_clean()
  end subroutine finish_homo_resistive_test






end module testmod_homo_resistive
