!
! MODULE: mod_equilibrium
!
!> @author
!> Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Module containing all equilibrium arrays.
!
module mod_equilibrium
  use mod_units
  use mod_types, only: density_type, temperature_type, bfield_type, velocity_type, &
                       gravity_type, resistivity_type, cooling_type, conduction_type
  use mod_global_variables, only: dp, gauss_gridpts, x_start, x_end, &
                                  flow, resistivity, external_gravity, radiative_cooling, &
                                  thermal_conduction, geometry, use_defaults, cgs_units
  use mod_physical_constants, only: dpi
  use mod_grid, only: initialise_grid, grid_gauss
  use mod_equilibrium_params, only: k2, k3
  use mod_logging, only: log_message, dp_fmt, char_log
  implicit none

  private

  procedure(), pointer :: set_equilibrium_values => null()

  interface
    module subroutine adiabatic_homo_eq; end subroutine
    module subroutine resistive_homo_eq; end subroutine
    module subroutine gravito_acoustic_eq; end subroutine
    module subroutine gravito_mhd_eq; end subroutine
    module subroutine resistive_tearing_modes_eq; end subroutine
    module subroutine resistive_tearing_modes_flow_eq; end subroutine
    module subroutine flow_driven_instabilities_eq; end subroutine
    module subroutine suydam_cluster_eq; end subroutine
    module subroutine kh_instability_eq; end subroutine
    module subroutine rotating_plasma_cyl_eq; end subroutine
    module subroutine kh_cd_instability_eq; end subroutine
    module subroutine internal_kink_eq; end subroutine
    module subroutine rotating_theta_pinch_eq; end subroutine
    module subroutine ideal_quasimodes_eq; end subroutine
    module subroutine uniform_thermal_cond_eq; end subroutine
    module subroutine nonuniform_thermal_cond_eq; end subroutine
    module subroutine magneto_rotational_eq; end subroutine
    module subroutine interface_modes_eq; end subroutine
    module subroutine discrete_alfven_eq; end subroutine
    module subroutine interchange_modes_eq; end subroutine
    module subroutine constant_current_eq; end subroutine
    module subroutine resonant_absorption_eq; end subroutine
    module subroutine magnetothermal_instability_eq; end subroutine
    module subroutine photospheric_flux_tube_eq; end subroutine
    module subroutine coronal_flux_tube_eq; end subroutine
    module subroutine gold_hoyle_eq; end subroutine

    module subroutine user_defined_eq; end subroutine

    module subroutine beta0_test_eq; end subroutine
    module subroutine hydro_test_eq; end subroutine
  end interface

  type (density_type)     :: rho_field
  type (temperature_type) :: T_field
  type (bfield_type)      :: B_field
  type (velocity_type)    :: v_field
  type (gravity_type)     :: grav_field
  type (resistivity_type) :: eta_field
  type (cooling_type)     :: rc_field
  type (conduction_type)  :: kappa_field

  public :: rho_field
  public :: T_field
  public :: B_field
  public :: v_field
  public :: grav_field
  public :: eta_field
  public :: rc_field
  public :: kappa_field

  public :: initialise_equilibrium
  public :: set_equilibrium
  public :: allow_geometry_override
  public :: equilibrium_clean


contains

  !> Initialises the equilibrium by allocating all equilibrium arrays and
  !! setting them to zero.
  subroutine initialise_equilibrium()
    use mod_types, only: initialise_type

    ! allocate and initialise everything
    call initialise_type(rho_field)
    call initialise_type(T_field)
    call initialise_type(B_field)
    call initialise_type(v_field)
    call initialise_type(grav_field)
    call initialise_type(eta_field)
    call initialise_type(rc_field)
    call initialise_type(kappa_field)

  end subroutine initialise_equilibrium


  subroutine set_equilibrium()
    use mod_check_values, only: check_negative_array, check_equilibrium_conditions, &
                                check_nan_values
    use mod_resistivity, only: set_resistivity_values
    use mod_radiative_cooling, only: initialise_radiative_cooling, set_radiative_cooling_values
    use mod_thermal_conduction, only: set_conduction_values

    ! Set equilibrium submodule to use
    call set_equilibrium_pointer()

    ! Call submodule
    call set_equilibrium_values()

    ! Set normalisations if needed
    call check_if_normalisations_set()

    ! Check equilibrium values, should be done before adding physics
    call check_negative_array(rho_field % rho0, 'density')
    call check_negative_array(T_field % T0, 'temperature')
    call check_equilibrium_conditions(rho_field, T_field, B_field, v_field, &
                                      grav_field, rc_field, kappa_field)
    call check_nan_values(rho_field)
    call check_nan_values(T_field)
    call check_nan_values(B_field)
    call check_nan_values(v_field)
    call check_nan_values(grav_field)

    ! Setup additional physics
    if (resistivity) then
      call set_resistivity_values(T_field, eta_field)
    end if
    if (radiative_cooling) then
      call initialise_radiative_cooling()
      call set_radiative_cooling_values(rho_field, T_field, rc_field)
    end if
    if (thermal_conduction) then
      call set_conduction_values(rho_field, T_field, B_field, kappa_field)
    end if

  end subroutine set_equilibrium


  subroutine set_equilibrium_pointer()
    use mod_global_variables, only: equilibrium_type

    select case(equilibrium_type)
    case('adiabatic_homo')
      set_equilibrium_values => adiabatic_homo_eq
    case('resistive_homo')
      set_equilibrium_values => resistive_homo_eq
    case("gravito_acoustic")
      set_equilibrium_values => gravito_acoustic_eq
    case("gravito_mhd")
      set_equilibrium_values => gravito_mhd_eq
    case("resistive_tearing")
      set_equilibrium_values => resistive_tearing_modes_eq
    case("resistive_tearing_flow")
      set_equilibrium_values => resistive_tearing_modes_flow_eq
    case("flow_driven_instabilities")
      set_equilibrium_values => flow_driven_instabilities_eq
    case("suydam_cluster")
      set_equilibrium_values => suydam_cluster_eq
    case("kelvin_helmholtz")
      set_equilibrium_values => kh_instability_eq
    case("rotating_plasma_cylinder")
      set_equilibrium_values => rotating_plasma_cyl_eq
    case("kelvin_helmholtz_cd")
      set_equilibrium_values => kh_cd_instability_eq
    case("internal_kink")
      set_equilibrium_values => internal_kink_eq
    case("rotating_theta_pinch")
      set_equilibrium_values => rotating_theta_pinch_eq
    case("ideal_quasimodes")
      set_equilibrium_values => ideal_quasimodes_eq
    case("uniform_conduction")
      set_equilibrium_values => uniform_thermal_cond_eq
    case("nonuniform_conduction")
      set_equilibrium_values => nonuniform_thermal_cond_eq
    case("magneto_rotational")
      set_equilibrium_values => magneto_rotational_eq
    case("interface_modes")
      set_equilibrium_values => interface_modes_eq
    case("discrete_alfven")
      set_equilibrium_values => discrete_alfven_eq
    case("interchange_modes")
      set_equilibrium_values => interchange_modes_eq
    case("constant_current_tokamak")
      set_equilibrium_values => constant_current_eq
    case("resonant_absorption")
      set_equilibrium_values => resonant_absorption_eq
    case("magnetothermal_instabilities")
      set_equilibrium_values => magnetothermal_instability_eq
    case("photospheric_flux_tube")
      set_equilibrium_values => photospheric_flux_tube_eq
    case("coronal_flux_tube")
      set_equilibrium_values => coronal_flux_tube_eq
    case("gold_hoyle")
      set_equilibrium_values => gold_hoyle_eq

    ! User defined
  case("user_defined")
      set_equilibrium_values => user_defined_eq

    ! Tests
  case("beta0_test")
      set_equilibrium_values => beta0_test_eq
    case("hydro_test")
      set_equilibrium_values => hydro_test_eq

    case default
      call log_message("equilibrium not recognised: " // trim(equilibrium_type), level='error')
    end select
  end subroutine set_equilibrium_pointer


  subroutine allow_geometry_override(default_geometry, default_x_start, default_x_end)
    use, intrinsic  :: ieee_arithmetic, only: ieee_is_nan
    use mod_global_variables, only: dp_LIMIT

    character(*), intent(in)  :: default_geometry
    real(dp), intent(in)      :: default_x_start, default_x_end

    if (geometry /= "" .and. geometry /= default_geometry) then
      call log_message("overriding default geometry with " // trim(geometry), level='warning')
    else
      geometry = default_geometry
    end if

    if ( (.not. ieee_is_nan(x_start)) .and. abs(x_start - default_x_start) >= dp_LIMIT ) then
      write(char_log, dp_fmt) x_start
      call log_message("overriding x_start with " // trim(char_log), level='warning')
    else
      x_start = default_x_start
    end if

    if ( (.not. ieee_is_nan(x_end)) .and. abs(x_end - default_x_end) >= dp_LIMIT ) then
      write(char_log, dp_fmt) x_end
      call log_message("overriding x_end with " // trim(char_log), level='warning')
    else
      x_end = default_x_end
    end if
  end subroutine allow_geometry_override


  !> Cleaning routine, deallocates all arrays in this module.
  subroutine equilibrium_clean()
    use mod_types, only: deallocate_type

    call deallocate_type(rho_field)
    call deallocate_type(T_field)
    call deallocate_type(B_field)
    call deallocate_type(v_field)
    call deallocate_type(grav_field)
    call deallocate_type(eta_field)
    call deallocate_type(rc_field)
    call deallocate_type(kappa_field)
  end subroutine equilibrium_clean

end module mod_equilibrium
