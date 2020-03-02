!
! PROGRAM: legolas
!
!> @author: Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Main program. Finite element code to calculate eigenvalues
!! and eigenvectors of the complete non-adiabatic MHD spectrum.
!! Included physics: flow, radiative cooling, thermal conduction, resistivity
program legolas
  use mod_global_variables, only: dp, str_len, show_results, run_silent
  use mod_matrix_creation, only: create_matrices
  use mod_solvers, only: solve_QR
  implicit none

  !> A matrix in eigenvalue problem wBX = AX
  complex(dp), allocatable  :: matrix_A(:, :)
  !> B matrix in eigenvalue problem wBX = AX
  real(dp), allocatable     :: matrix_B(:, :)
  !> Solutions to the eigenvalue problem
  complex(dp), allocatable  :: omega(:)
  !> Right eigenvectors
  complex(dp), allocatable  :: eigenvecs_right(:, :)
  !> Left eigenvectors
  complex(dp), allocatable  :: eigenvecs_left(:, :)
  !> Name of the configuration file
  character(str_len)        :: config_file

  ! allocate variables, initialise grid and physics, set equilibrium
  call initialisation()

  ! create matrices A and B
  call create_matrices(matrix_B, matrix_A)

  ! solve eigenvalue problem
  if (.not. run_silent) then
    write(*, *) "Solving eigenvalue problem..."
  end if
  call solve_QR(matrix_A, matrix_B, omega, eigenvecs_left, eigenvecs_right)

  ! write spectrum, eigenvectors, matrices etc. to file
  call finalise_results()

  ! deallocate everything
  call cleanup()

  ! call python script and pass configuration file
  if (show_results) then
    write(*, *) ""
    write(*, *) "Plotting results..."
    call execute_command_line("python3 python/legolas_analyser.py " // &
                              "-i " // trim(config_file))
  end if

contains

  !> Initialises the grid and equilibrium configuration.
  subroutine initialisation()
    use mod_global_variables, only: matrix_gridpts
    use mod_input, only: read_parfile, get_parfile
    use mod_equilibrium, only: initialise_equilibrium, set_equilibrium
    use mod_eigenfunctions, only: initialise_eigenfunctions
    use mod_output, only: startup_info_toconsole

    character(len=str_len)  :: parfile

    ! print empty line
    write(*, *) ""

    call get_parfile(parfile)
    call read_parfile(parfile)

    ! Allocate matrices
    allocate(matrix_A(matrix_gridpts, matrix_gridpts))
    allocate(matrix_B(matrix_gridpts, matrix_gridpts))
    allocate(omega(matrix_gridpts))
    allocate(eigenvecs_right(matrix_gridpts, matrix_gridpts))
    allocate(eigenvecs_left(matrix_gridpts, matrix_gridpts))

    ! Initialise equilibrium
    call initialise_equilibrium()

    ! Initialise eigenfunction arrays
    call initialise_eigenfunctions()

    ! set the equilibrium configuration
    call set_equilibrium()

    ! print configuration to console
    if (.not. run_silent) then
      call startup_info_toconsole()
    end if

  end subroutine initialisation


  !> Saves the solutions
  subroutine finalise_results()
    use mod_global_variables, only: savename_eigenvalues, write_matrices, savename_matrix, write_eigenvectors, &
                                    savename_eigenvectors, write_eigenfunctions, write_equilibrium, &
                                    savename_equil, savename_config
    use mod_output, only: eigenvalues_tofile, matrices_tofile, eigenvectors_tofile, &
                          configuration_tofile, equilibrium_tofile
    use mod_eigenfunctions, only: calculate_eigenfunctions
    use mod_equilibrium, only: rho_field, T_field, B_field, v_field
    use mod_grid, only: grid_gauss

    call eigenvalues_tofile(omega, savename_eigenvalues)

    if (write_matrices) then
      call matrices_tofile(matrix_A, matrix_B, savename_matrix)
    end if

    if (write_eigenvectors) then
      call eigenvectors_tofile(eigenvecs_left, eigenvecs_right, savename_eigenvectors)
    end if

    if (write_eigenfunctions) then
      call calculate_eigenfunctions(eigenvecs_right)
    end if

    if (write_equilibrium) then
      call equilibrium_tofile(grid_gauss, rho_field, T_field, B_field, v_field, savename_equil)
    end if

    ! save running configuration
    call configuration_tofile(savename_config, config_file)
  end subroutine finalise_results


  !> Performs cleanup, deallocates variables.
  subroutine cleanup()
    use mod_global_variables, only: radiative_cooling
    use mod_grid, only: grid_clean
    use mod_equilibrium, only: equilibrium_clean
    use mod_radiative_cooling, only: radiative_cooling_clean
    use mod_eigenfunctions, only: eigenfunctions_clean

    deallocate(matrix_A)
    deallocate(matrix_B)
    deallocate(omega)
    deallocate(eigenvecs_left)
    deallocate(eigenvecs_right)

    call grid_clean()
    call equilibrium_clean()

    if (radiative_cooling) then
      call radiative_cooling_clean()
    end if
    call eigenfunctions_clean()

  end subroutine cleanup

end program legolas
