module mod_io
  use mod_global_variables
  use mod_physical_constants
  implicit none


  public

  integer, parameter  :: w_output = 10
  integer, parameter  :: config   = 20

contains

  subroutine save_eigenvalues(omega)
    complex(dp), intent(in) :: omega(matrix_gridpts)
    integer                 :: i

    open (w_output, file='eigenvalues.txt', status='unknown')

    do i = 1, matrix_gridpts
      write(w_output, *) omega(i)
    end do

  end subroutine save_eigenvalues

  subroutine save_config()
    open(config, file='config.txt', status='unknown')

    write(config, *) "Equilibrium type   : ", equilibrium_type

    write(config, *) "Flow               : ", flow
    write(config, *) "Radiative cooling  : ", radiative_cooling
    if (radiative_cooling) then
      write(config, *) "Cooling curve      : ", cooling_curve
    end if
    write(config, *) "Thermal conduction : ", thermal_conduction
    write(config, *) "Resistivity        : ", resistivity
    write(config, *) "External gravity   : ", external_gravity
    if (external_gravity) then
      write(config, *) "Gravity strength   : ", gravity_type
    end if

    write(config, *) ""

    write(config, *) "Geometry           : ", geometry
    write(config, *) "Start              : ", x_start
    write(config, *) "End                : ", x_end
    write(config, *) "Gridpoints         : ", gridpts
    write(config, *) "Matrix gridpoints  : ", matrix_gridpts
    write(config, *) "Gaussian gridpoints: ", gauss_gridpts
    write(config, *) "Mesh accumulation  : ", mesh_accumulation

    write(config, *) ""

    write(config, *) "CGS units          : ", cgs_units
    write(config, *) "Unit length        : ", unit_length
    write(config, *) "Unit time          : ", unit_time
    write(config, *) "Unit density       : ", unit_density
    write(config, *) "Unit velocity      : ", unit_velocity
    write(config, *) "Unit temperature   : ", unit_temperature
    write(config, *) "Unit pressure      : ", unit_pressure
    write(config, *) "Unit magnetic field: ", unit_magneticfield
    write(config, *) "Unit numberdensity : ", unit_numberdensity
    write(config, *) "Unit luminosity    : ", unit_luminosity
    write(config, *) "Unit resistivity   : ", unit_resistivity
  end subroutine save_config

end module mod_io
