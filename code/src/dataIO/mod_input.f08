module mod_input
  use mod_global_variables
  use mod_physical_constants
  implicit none

contains

  subroutine read_parfile()
    real(dp)    :: mhd_gamma
    real(dp)    :: unit_length, unit_numberdensity, unit_temperature, &
                   unit_velocity
    integer     :: gridpoints

    namelist /gridlist/ geometry, x_start, x_end, gridpoints
    namelist /meshlist/ mesh_accumulation, ev_1, ev_2, sigma_1, sigma_2
    namelist /physicslist/ mhd_gamma, flow, radiative_cooling, ncool, &
                           cooling_curve, external_gravity, gravity_type, &
                           thermal_conduction, resistivity
    namelist /unitslist/ cgs_units, unit_length, unit_numberdensity, &
                         unit_temperature, unit_velocity
    namelist /equilibriumlist/ use_precoded, equilibrium_type, boundary_type
    namelist /savelist/ plot_when_finished, write_AB, write_eigenvectors, &
                        write_eigenfunctions

    !! Set defaults
    !> Gridlist defaults
    geometry = "Cartesian"          !< geometry of the problem
    x_start  = 0.0d0                !< start of the grid
    x_end    = 1.0d0                !< end of the grid
    gridpoints = 11

    !> Meshlist defaults
    mesh_accumulation = .false.
    ev_1 = 1.25d0                   !< expected value gaussian 1
    ev_2 = 1.25d0                   !< expected value gaussian 2
    sigma_1 = 1.0d0                 !< standard deviation gaussian 1
    sigma_2 = 2.0d0                 !< standard deviation gaussian 2

    !> Physicslist defaults
    mhd_gamma = 5.0d0 / 3.0d0       !< ratio of specific heats
    flow  = .false.                 !< use flow
    radiative_cooling = .false.     !< use radiative cooling
    ncool = 4000                    !< amount of points to interpolate curve
    cooling_curve = "JCcorona"      !< cooling curve to use
    external_gravity = .false.      !< use external gravity
    gravity_type = "solar"          !< strength of external gravity
    thermal_conduction = .false.    !< use thermal conduction
    resistivity = .false.           !< use resistivity

    !> Unitslist defaults
    cgs_units = .true.              !< use cgs units
    unit_length = 1.0d9
    unit_numberdensity = 1.0d9
    unit_temperature = 1.0d6
    unit_velocity = 0.0d0

    !> Equilibriumlist defaults
    use_precoded = .true.                       !< use precoded equilibrium
    equilibrium_type = "Adiabatic homogeneous"
    !equilibrium_type = "Resistive homogeneous"
    !equilibrium_type = "Suydam cluster modes"
    !equilibrium_type = "Kelvin-Helmholtz"
    boundary_type = 'wall'

    !> Savelist defaults
    plot_when_finished = .true.     !< plot spectrum when finished
    write_AB = .true.              !< write matrices A and B when finished
    write_eigenvectors = .true.    !< writes eigenvectors to file
    write_eigenfunctions = .true.   !< writes eigenfunctions to file


    call set_gridpts(gridpoints)
    call set_gamma(mhd_gamma)

    call set_unit_length(unit_length)
    call set_unit_numberdensity(unit_numberdensity)
    call set_unit_temperature(unit_temperature)
    call set_normalisations(cgs_units)

  end subroutine read_parfile

end module mod_input
