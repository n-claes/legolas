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
  use mod_types, only: density_type, temperature_type, bfield_type, velocity_type, &
                       gravity_type, resistivity_type, cooling_type, conduction_type
  use mod_global_variables, only: dp, gauss_gridpts, k2, k3, x_start, x_end, &
                                  flow, resistivity, external_gravity, radiative_cooling, &
                                  thermal_conduction, geometry
  use mod_physical_constants, only: dpi
  use mod_grid, only: initialise_grid, grid_gauss
  implicit none

  private

  abstract interface
    subroutine eq_void()
    end subroutine eq_void
  end interface

  procedure (eq_void), pointer :: set_equilibrium_values => null()

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
    module subroutine tokamak_cyl_eq; end subroutine
    module subroutine nonuniform_thermal_cond_eq; end subroutine
    module subroutine magneto_rotational_eq; end subroutine
    module subroutine interface_modes_eq; end subroutine
    module subroutine discrete_alfven_eq; end subroutine
    module subroutine interchange_modes_eq; end subroutine

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
  public :: equilibrium_clean


contains

  !> Initialises the equilibrium by allocating all equilibrium arrays and
  !! setting them to zero.
  subroutine initialise_equilibrium()
    use mod_types, only: initialise_type
    use mod_radiative_cooling, only: initialise_radiative_cooling

    ! allocate and initialise everything
    call initialise_type(rho_field)
    call initialise_type(T_field)
    call initialise_type(B_field)
    call initialise_type(v_field)
    call initialise_type(grav_field)
    call initialise_type(eta_field)
    call initialise_type(rc_field)
    call initialise_type(kappa_field)

    ! initialise radiative cooling curves
    if (radiative_cooling) then
      call initialise_radiative_cooling()
    end if

  end subroutine initialise_equilibrium


  subroutine set_equilibrium()
    use mod_check_values, only: check_negative_array, check_equilibrium_conditions
    use mod_resistivity, only: set_resistivity_values
    use mod_radiative_cooling, only: set_radiative_cooling_values
    use mod_thermal_conduction, only: set_conduction_values

    ! Set equilibrium submodule to use
    call set_equilibrium_pointer()
    ! Call submodule
    call set_equilibrium_values()

    ! Setup additional physics
    if (resistivity) then
      call set_resistivity_values(T_field, eta_field)
    end if
    if (radiative_cooling) then
      call set_radiative_cooling_values(rho_field, T_field, rc_field)
    end if
    if (thermal_conduction) then
      call set_conduction_values(rho_field, T_field, B_field, kappa_field)
    end if

    ! Check equilibrium values
    call check_negative_array(rho_field % rho0, 'density')
    call check_negative_array(T_field % T0, 'temperature')
    call check_equilibrium_conditions(rho_field, T_field, B_field, v_field, grav_field)

  end subroutine set_equilibrium


  subroutine set_equilibrium_pointer()
    use mod_global_variables, only: equilibrium_type

    select case(equilibrium_type)
    case('Adiabatic homogeneous')
      set_equilibrium_values => adiabatic_homo_eq
    case('Resistive homogeneous')
      set_equilibrium_values => resistive_homo_eq
    case("Gravito-acoustic waves")
      set_equilibrium_values => gravito_acoustic_eq
    case("Gravito-MHD waves")
      set_equilibrium_values => gravito_mhd_eq
    case("Resistive tearing modes")
      set_equilibrium_values => resistive_tearing_modes_eq
    case("Resistive tearing modes with flow")
      set_equilibrium_values => resistive_tearing_modes_flow_eq
    case("Flow driven instabilities")
      set_equilibrium_values => flow_driven_instabilities_eq
    case("Suydam cluster modes")
      set_equilibrium_values => suydam_cluster_eq
    case("Kelvin-Helmholtz")
      set_equilibrium_values => kh_instability_eq
    case("Rotating plasma cylinder")
      set_equilibrium_values => rotating_plasma_cyl_eq
    case("Kelvin-Helmholtz and current driven")
      set_equilibrium_values => kh_cd_instability_eq
    case("Internal kink modes")
      set_equilibrium_values => internal_kink_eq
    case("Rotating theta pinch")
      set_equilibrium_values => rotating_theta_pinch_eq
    case("Ideal quasimodes")
      set_equilibrium_values => ideal_quasimodes_eq
    case("Simple tokamak with resistivity")
      set_equilibrium_values => tokamak_cyl_eq
    case("Uniform with thermal conduction")
      set_equilibrium_values => uniform_thermal_cond_eq
    case("Non-uniform with thermal conduction")
      set_equilibrium_values => nonuniform_thermal_cond_eq
    case("Magneto-rotational instability")
      set_equilibrium_values => magneto_rotational_eq
    case("Interface modes")
      set_equilibrium_values => interface_modes_eq
    case("Non-adiabatic discrete Alfven")
      set_equilibrium_values => discrete_alfven_eq
    case("Interchange modes")
      set_equilibrium_values => interchange_modes_eq

    ! User defined
    case("User defined equilibrium")
      set_equilibrium_values => user_defined_eq

    ! Tests
    case("Beta=0 test")
      set_equilibrium_values => beta0_test_eq
    case("Hydrodynamics test")
      set_equilibrium_values => hydro_test_eq

    case default
      write(*, *) "Equilibrium not recognised."
      write(*, *) "currently set on: ", equilibrium_type
      stop
    end select
  end subroutine set_equilibrium_pointer


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
