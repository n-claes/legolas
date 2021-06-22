! =============================================================================
!> Parent module governing all equilibrium types and submodules.
!! This module contains all equilibrium types and the initial
!! declarations of the module subroutines. Every equilibrium submodule
!! extends this module, implementing one of the module subroutines declared here.
!! All "main" equilibrium configurations are set in the submodules. The ones that
!! depend on "main" arrays, like radiative cooling, are set here through calls to their
!! respective modules.
!! @note    All use statements specified here at the main module scope
!!          are automatically accessible in every submodule that extends this one.
module mod_equilibrium
  use mod_units
  use mod_types
  use mod_global_variables, only: dp, gauss_gridpts, x_start, x_end, &
    flow, resistivity, external_gravity, radiative_cooling, &
    thermal_conduction, viscosity, hall_mhd, geometry, use_defaults, cgs_units
  use mod_physical_constants, only: dpi
  use mod_grid, only: initialise_grid, grid_gauss
  use mod_equilibrium_params, only: k2, k3
  use mod_logging, only: log_message, str
  implicit none

  private

  !> pointer for the submodule, initialised to null
  procedure(), pointer :: set_equilibrium_values => null()

  !> interface to the different equilibrium submodules
  interface
    module subroutine adiabatic_homo_eq; end subroutine
    module subroutine constant_current_eq; end subroutine
    module subroutine coronal_flux_tube_eq; end subroutine
    module subroutine discrete_alfven_eq; end subroutine
    module subroutine flow_driven_instabilities_eq; end subroutine
    module subroutine gold_hoyle_eq; end subroutine
    module subroutine gravito_acoustic_eq; end subroutine
    module subroutine gravito_mhd_eq; end subroutine
    module subroutine interchange_modes_eq; end subroutine
    module subroutine internal_kink_eq; end subroutine
    module subroutine isothermal_atmosphere_eq; end subroutine
    module subroutine KHI_eq; end subroutine
    module subroutine kh_cd_instability_eq; end subroutine
    module subroutine magnetothermal_instability_eq; end subroutine
    module subroutine MRI_accretion_eq; end subroutine
    module subroutine photospheric_flux_tube_eq; end subroutine
    module subroutine resistive_homo_eq; end subroutine
    module subroutine resistive_tearing_modes_eq; end subroutine
    module subroutine resistive_tearing_modes_flow_eq; end subroutine
    module subroutine resonant_absorption_eq; end subroutine
    module subroutine rotating_plasma_cyl_eq; end subroutine
    module subroutine RTI_eq; end subroutine
    module subroutine RTI_KHI_eq; end subroutine
    module subroutine RTI_theta_pinch_eq; end subroutine
    module subroutine suydam_cluster_eq; end subroutine
    module subroutine viscoresistive_tube_eq; end subroutine
    module subroutine couette_flow_eq; end subroutine
    module subroutine taylor_couette_eq; end subroutine
    module subroutine harris_sheet_eq; end subroutine
    module subroutine tc_pinch_eq; end subroutine
    module subroutine tc_hmri_eq; end subroutine
    module subroutine user_defined_eq; end subroutine
  end interface

  !> type containing all density-related equilibrium variables
  type (density_type)     :: rho_field
  !> type containing all temperature-related equilibrium variables
  type (temperature_type) :: T_field
  !> type containing all magnetic field-related equilibrium variables
  type (bfield_type)      :: B_field
  !> type containing all velocity-related equilibrium variables
  type (velocity_type)    :: v_field
  !> type containing all gravity-related equilibrium variables
  type (gravity_type)     :: grav_field
  !> type containing all resistivity-related equilibrium variables
  type (resistivity_type) :: eta_field
  !> type containig all radiative cooling-related equilibrium variables
  type (cooling_type)     :: rc_field
  !> type containing all thermal conduction-related equilibrium variables
  type (conduction_type)  :: kappa_field
  !> type containing all Hall related variables
  type (hall_type)        :: hall_field

  public :: rho_field
  public :: T_field
  public :: B_field
  public :: v_field
  public :: grav_field
  public :: eta_field
  public :: rc_field
  public :: kappa_field
  public :: hall_field

  public :: initialise_equilibrium
  public :: set_equilibrium
  public :: allow_geometry_override
  public :: equilibrium_clean

contains


  !> Initialises the equilibrium types by calling the corresponding
  !! subroutine, which allocates all necessary attributes.
  subroutine initialise_equilibrium()
    call initialise_type(rho_field)
    call initialise_type(T_field)
    call initialise_type(B_field)
    call initialise_type(v_field)
    call initialise_type(grav_field)
    call initialise_type(eta_field)
    call initialise_type(rc_field)
    call initialise_type(kappa_field)
    call initialise_type(hall_field)
  end subroutine initialise_equilibrium


  !> Calls the routine to set the equilibrium pointer, then calls the correct
  !! submodule. Performs some sanity checks (negative values, NaNs etc.) when
  !! the equilibrium is set, then calls additional physics modules if needed.
  !! @warning Throws appropriate errors if the equilibrium configuration is
  !!          not balanced, contains NaN or if density/temperature contains
  !!          negative values.
  subroutine set_equilibrium()
    use mod_global_variables, only: coaxial, dp_LIMIT
    use mod_inspections, only: perform_NaN_and_negative_checks, perform_sanity_checks
    use mod_resistivity, only: set_resistivity_values
    use mod_radiative_cooling, only: initialise_radiative_cooling, &
      set_radiative_cooling_values
    use mod_thermal_conduction, only: set_conduction_values
    use mod_hall, only: set_hall_factors

    ! Set equilibrium submodule to use
    call set_equilibrium_pointer()
    ! Call submodule
    call set_equilibrium_values()
    ! Set normalisations if needed
    call check_if_normalisations_set()

    ! Check x_start if coaxial is true
    if (coaxial .and. x_start <= dp_LIMIT) then
      call log_message( &
        "x_start must be > 0 to introduce an inner wall boundary", level="error" &
      )
      return
    end if

    ! Do initial checks for NaN and negative density/temperature
    call perform_NaN_and_negative_checks( &
      rho_field, T_field, B_field, v_field, grav_field &
    )

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
    if (hall_mhd) then
      call set_hall_factors(hall_field)
    end if

    ! Do final sanity checks on values
    call perform_sanity_checks( &
      rho_field, T_field, B_field, v_field, grav_field, rc_field, kappa_field &
    )
  end subroutine set_equilibrium


  !> Selects the submodule based on the specified equilibrium
  !! in the parfile. Works on a case-select basis.
  !! @warning   Throws an error if the equilibrium type is not recognised.
  subroutine set_equilibrium_pointer()
    use mod_global_variables, only: equilibrium_type

    select case(equilibrium_type)
    case("adiabatic_homo")
      set_equilibrium_values => adiabatic_homo_eq
    case("constant_current_tokamak")
      set_equilibrium_values => constant_current_eq
    case("coronal_flux_tube")
      set_equilibrium_values => coronal_flux_tube_eq
    case("discrete_alfven")
      set_equilibrium_values => discrete_alfven_eq
    case("gold_hoyle")
      set_equilibrium_values => gold_hoyle_eq
    case("gravito_acoustic")
      set_equilibrium_values => gravito_acoustic_eq
    case("gravito_mhd")
      set_equilibrium_values => gravito_mhd_eq
    case("interchange_modes")
      set_equilibrium_values => interchange_modes_eq
    case("internal_kink")
      set_equilibrium_values => internal_kink_eq
    case("isothermal_atmosphere")
      set_equilibrium_values => isothermal_atmosphere_eq
    case("kelvin_helmholtz")
      set_equilibrium_values => KHI_eq
    case("kelvin_helmholtz_cd")
      set_equilibrium_values => kh_cd_instability_eq
    case("magnetothermal_instabilities")
      set_equilibrium_values => magnetothermal_instability_eq
    case("MRI_accretion")
      set_equilibrium_values => MRI_accretion_eq
    case("photospheric_flux_tube")
      set_equilibrium_values => photospheric_flux_tube_eq
    case("rayleigh_taylor")
      set_equilibrium_values => RTI_eq
    case("resistive_homo")
      set_equilibrium_values => resistive_homo_eq
    case("resistive_tearing")
      set_equilibrium_values => resistive_tearing_modes_eq
    case("resistive_tearing_flow")
      set_equilibrium_values => resistive_tearing_modes_flow_eq
    case("resonant_absorption")
      set_equilibrium_values => resonant_absorption_eq
    case("rotating_plasma_cylinder")
      set_equilibrium_values => rotating_plasma_cyl_eq
    case("RTI_KHI")
      set_equilibrium_values => RTI_KHI_eq
    case("RTI_theta_pinch")
      set_equilibrium_values => RTI_theta_pinch_eq
    case("suydam_cluster")
      set_equilibrium_values => suydam_cluster_eq
    case("viscoresistive_tube")
      set_equilibrium_values => viscoresistive_tube_eq
    case("couette_flow")
      set_equilibrium_values => couette_flow_eq
    case("taylor_couette")
      set_equilibrium_values => taylor_couette_eq
    case("harris_sheet")
      set_equilibrium_values => harris_sheet_eq
    case("tc_pinch")
      set_equilibrium_values => tc_pinch_eq
    case("tc_HMRI")
      set_equilibrium_values => tc_hmri_eq
    case("user_defined")
      set_equilibrium_values => user_defined_eq
    case default
      call log_message( &
        "equilibrium not recognised: " // trim(equilibrium_type), level="error" &
      )
    end select
  end subroutine set_equilibrium_pointer


  !> Allows overriding geometry and grid-related parameters.
  !! Sets default values for the geometry and grid start/end. If this subroutine is
  !! used to set geometry/grid values in the submodule it becomes possible to override
  !! them through the parfile. Warnings will always be printed if this happens.
  !! @note  If values specified in the parfile are equal to the default values,
  !!        nothing happens. @endnote
  !! @warning A warning is thrown if:
  !!
  !! - the geometry specified in the submodule is overridden.
  !! - the starting value of the grid specified in the submodule is overridden.
  !! - the end value of the grid specified in the submodule is overridden. @endwarning
  subroutine allow_geometry_override(default_geometry, default_x_start, default_x_end)
    use mod_check_values, only: is_NaN
    use mod_global_variables, only: dp_LIMIT

    !> default geometry to set
    character(*), intent(in), optional  :: default_geometry
    !> default start of the grid
    real(dp), intent(in), optional  :: default_x_start
    !> default end of the grid
    real(dp), intent(in), optional  :: default_x_end

    if (present(default_geometry)) then
      if (geometry /= "" .and. geometry /= default_geometry) then
        call log_message( &
          "overriding default geometry with " // trim(geometry), level="warning" &
        )
      else
        geometry = default_geometry
      end if
    end if

    if (present(default_x_start)) then
      if ( &
        (.not. is_NaN(x_start)) .and. abs(x_start - default_x_start) >= dp_LIMIT &
      ) then
        call log_message("overriding x_start with " // str(x_start), level="warning")
      else
        x_start = default_x_start
      end if
    end if

    if (present(default_x_end)) then
      if ( &
        (.not. is_NaN(x_end)) .and. abs(x_end - default_x_end) >= dp_LIMIT &
      ) then
        call log_message("overriding x_end with " // str(x_end), level="warning")
      else
        x_end = default_x_end
      end if
    end if
  end subroutine allow_geometry_override


  !> Cleaning routine, deallocates the equilibrium types.
  subroutine equilibrium_clean()
    call deallocate_type(rho_field)
    call deallocate_type(T_field)
    call deallocate_type(B_field)
    call deallocate_type(v_field)
    call deallocate_type(grav_field)
    call deallocate_type(eta_field)
    call deallocate_type(rc_field)
    call deallocate_type(kappa_field)
    call deallocate_type(hall_field)
  end subroutine equilibrium_clean

end module mod_equilibrium
