module testmod_homogeneous
  use mod_global_variables
  use mod_grid
  use mod_equilibrium
  use mod_equilibrium_derivatives
  use mod_setup_matrix_a
  use mod_setup_matrix_b

  implicit none


  real(dp), allocatable        :: mat_B(:, :)
  complex(dp), allocatable     :: mat_A(:, :)
  complex(dp), allocatable     :: omega(:)

  integer, parameter           :: gridpts_homo = 11
  real(dp), parameter          :: eps = 1.0d-14


contains

  subroutine init_homogeneous()
    use_precoded = .true.
    equilibrium_type = "Adiabatic homogeneous"
    !! Set global variables
    flow = .false.
    radiative_cooling = .false.
    thermal_conduction = .false.
    resistivity = .false.
    external_gravity = .false.

    !! Set grid points
    x_start = 0.0d0
    x_end   = 1.0d0
    call set_gridpts(gridpts_homo)
    call set_gamma(5.0d0 / 3.0d0)
    call initialise_grid()
    call set_grid_gauss()

    call initialise_equilibrium()
    call initialise_equilibrium_derivatives()
    call set_equilibrium()

    allocate(mat_B(matrix_gridpts, matrix_gridpts))
    allocate(mat_A(matrix_gridpts, matrix_gridpts))
    allocate(omega(matrix_gridpts))
  end subroutine init_homogeneous

  subroutine run_homogeneous_test()
    use mod_solvers
    use mod_physical_constants
    use mod_io

    complex(dp)   :: mat_A_theory(8, 8)
    real(dp)      :: mat_B_theory(8, 8)

    complex(dp)   :: omega_theory(8)
    complex(dp)   :: ic0
    real(dp)      :: k1, k2, k3, va, vs
    real(dp)      :: P0, T0, B03, rho0
    integer       :: n, it_end, i

    real(dp)      :: w_real, w_imag
    logical       :: append

    ic0 = (1.0d0, 0.0d0)

    !! Solve using the code itself
    call construct_B(mat_B)
    call construct_A(mat_A)
    call solve_QR(mat_A, mat_B, omega)

    !! Rely on analytical solution and quantization of wave number,
    !! kx will be quantized according to kx = n*pi / L
    mat_A_theory = 0.0d0
    mat_B_theory = 0.0d0

    ! Set variables (use 1 for index but does not matter as it is homogeneous)
    rho0 = rho0_eq(1)
    T0   = T0_eq(1)
    P0   = rho0 * T0
    B03  = B03_eq(1)

    ! construct matrix B
    mat_B_theory(1, 1) = gamma
    mat_B_theory(2, 2) = 1.0d0
    mat_B_theory(3, 3) = 1.0d0
    mat_B_theory(4, 4) = 1.0d0
    mat_B_theory(5, 5) = gamma / gamma_1
    mat_B_theory(6, 6) = 1.0d0
    mat_B_theory(7, 7) = 1.0d0
    mat_B_theory(8, 8) = 1.0d0

    ! sound speed   sqrt(gamma * P0 / rho0)
    vs = sqrt(gamma * P0 / rho0)
    ! Alfven speed  B0z / sqrt(rho0)    rho0 = 1, B0z = 1
    va = B03 / sqrt(rho0)

    k2 = 0.0d0
    k3 = dpi

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
      mat_A_theory(1, 2) = k1 * vs
      mat_A_theory(1, 4) = k3 * vs
      mat_A_theory(2, 1) = k1 * vs
      mat_A_theory(2, 5) = k1 * vs
      mat_A_theory(2, 6) = -k3 * va
      mat_A_theory(2, 8) = k1 * va
      mat_A_theory(3, 7) = -k3 * va
      mat_A_theory(4, 1) = k3 * vs
      mat_A_theory(4, 5) = k3 * vs
      mat_A_theory(5, 2) = k1 * vs
      mat_A_theory(5, 4) = k3 * vs
      mat_A_theory(6, 2) = -k3 * va
      mat_A_theory(7, 3) = -k3 * va
      mat_A_theory(8, 2) = k1 * va

      ! make it complex
      mat_A_theory = mat_A_theory * (1.0d0, 0.0d0)
      call solve_QR(mat_A_theory, mat_B_theory, omega_theory)

      ! For spectrum plotting, force very small numbers to 0
      do i = 1, 8
        w_real = real(omega_theory(i))
        w_imag = aimag(omega_theory(i))
        if (abs(w_real) < eps) then
          w_real = 0.0d0
        end if
        if (abs(w_imag) < eps) then
          w_imag = 0.0d0
        end if
        omega_theory(i) = cmplx(w_real, w_imag, kind=dp)
      end do

      !after writing once, start appending to file
      if (n == 2) then
        append = .true.
      end if

      call save_eigenvalues(omega_theory, "omega_test_homo", append)

    end do

  end subroutine run_homogeneous_test

  subroutine finish_homogeneous()
    deallocate(mat_B)
    deallocate(mat_A)
    deallocate(omega)
    write(*, *) "TEST COMPLETED"
    call grid_clean()
    call equilibrium_clean()
    call equilibrium_derivatives_clean()
  end subroutine finish_homogeneous






end module testmod_homogeneous
