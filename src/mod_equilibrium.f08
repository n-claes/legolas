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
  use mod_types
  use mod_global_variables, only: dp
  use mod_physical_constants, only: dpi
  use mod_grid, only: initialise_grid, grid_gauss
  use mod_equilibrium_params, only: k2, k3
  use mod_logging, only: logger, str, exp_fmt
  use mod_settings, only: settings_t
  use mod_background, only: background_t
  implicit none

  private

  !> pointer for the submodule, initialised to null
  procedure(), pointer :: set_equilibrium_values => null()

  !> interface to the different equilibrium submodules
  interface
    module subroutine adiabatic_homo_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine adiabatic_homo_eq
    module subroutine constant_current_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine constant_current_eq
    module subroutine coronal_flux_tube_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine coronal_flux_tube_eq
    module subroutine discrete_alfven_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine discrete_alfven_eq
    module subroutine flow_driven_instabilities_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine flow_driven_instabilities_eq
    module subroutine gold_hoyle_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine gold_hoyle_eq
    module subroutine gravito_acoustic_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine gravito_acoustic_eq
    module subroutine gravito_mhd_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine gravito_mhd_eq
    module subroutine interchange_modes_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine interchange_modes_eq
    module subroutine internal_kink_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine internal_kink_eq
    module subroutine isothermal_atmosphere_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine isothermal_atmosphere_eq
    module subroutine KHI_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine
    module subroutine kh_cd_instability_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine kh_cd_instability_eq
    module subroutine magnetothermal_instability_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine magnetothermal_instability_eq
    module subroutine MRI_accretion_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine MRI_accretion_eq
    module subroutine photospheric_flux_tube_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine photospheric_flux_tube_eq
    module subroutine resistive_homo_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine resistive_homo_eq
    module subroutine resistive_tearing_modes_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine resistive_tearing_modes_eq
    module subroutine resistive_tearing_modes_flow_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine resistive_tearing_modes_flow_eq
    module subroutine resonant_absorption_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine resonant_absorption_eq
    module subroutine rotating_plasma_cyl_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine rotating_plasma_cyl_eq
    module subroutine RTI_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine RTI_eq
    module subroutine RTI_KHI_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine RTI_KHI_eq
    module subroutine RTI_theta_pinch_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine RTI_theta_pinch_eq
    module subroutine suydam_cluster_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine suydam_cluster_eq
    module subroutine couette_flow_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine couette_flow_eq
    module subroutine taylor_couette_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine taylor_couette_eq
    module subroutine harris_sheet_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine harris_sheet_eq
    module subroutine tc_pinch_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine tc_pinch_eq
    module subroutine user_defined_eq(settings, background)
      type(settings_t), intent(inout) :: settings
      type(background_t), intent(inout) :: background
    end subroutine user_defined_eq
  end interface

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

  public :: grav_field
  public :: eta_field
  public :: rc_field
  public :: kappa_field
  public :: hall_field

  public :: initialise_equilibrium
  public :: set_equilibrium
  public :: equilibrium_clean

contains


  !> Initialises the equilibrium types by calling the corresponding
  !! subroutine, which allocates all necessary attributes.
  subroutine initialise_equilibrium(settings)
    type(settings_t), intent(inout) :: settings

    call initialise_type(settings, grav_field)
    call initialise_type(settings, eta_field)
    call initialise_type(settings, rc_field)
    call initialise_type(settings, kappa_field)
    call initialise_type(settings, hall_field)
  end subroutine initialise_equilibrium


  !> Calls the routine to set the equilibrium pointer, then calls the correct
  !! submodule. Performs some sanity checks (negative values, NaNs etc.) when
  !! the equilibrium is set, then calls additional physics modules if needed.
  !! @warning Throws appropriate errors if the equilibrium configuration is
  !!          not balanced, contains NaN or if density/temperature contains
  !!          negative values.
  subroutine set_equilibrium(settings, background)
    use mod_global_variables, only: dp_LIMIT
    use mod_inspections, only: perform_NaN_and_negative_checks, perform_sanity_checks
    use mod_resistivity, only: set_resistivity_values
    use mod_radiative_cooling, only: initialise_radiative_cooling, &
      set_radiative_cooling_values
    use mod_thermal_conduction, only: set_conduction_values
    use mod_hall, only: set_hall_factors
    type(settings_t), intent(inout) :: settings
    type(background_t), intent(inout) :: background

    ! Set equilibrium submodule to use
    call set_equilibrium_pointer(settings)
    ! Call submodule
    call set_equilibrium_values(settings, background)

    ! Check x_start if coaxial is true
    if (settings%grid%coaxial .and. settings%grid%get_grid_start() <= dp_LIMIT) then
      call logger%error("x_start must be > 0 to introduce an inner wall boundary")
      return
    end if

    ! Do initial checks for NaN and negative density/temperature
    call perform_NaN_and_negative_checks(background, grav_field)

    ! Setup additional physics
    if (settings%physics%resistivity%is_enabled()) then
      call set_resistivity_values(settings, background, eta_field)
    end if
    if (settings%physics%cooling%is_enabled()) then
      call initialise_radiative_cooling(settings)
      call set_radiative_cooling_values(settings, background, rc_field)
    end if
    if (settings%physics%conduction%is_enabled()) then
      call set_conduction_values(settings, background, kappa_field)
    end if
    if (settings%physics%hall%is_enabled()) then
      call set_hall_factors(settings, hall_field)
    end if

    ! Do final sanity checks on values
    call perform_sanity_checks(settings, background, grav_field, rc_field, kappa_field)
  end subroutine set_equilibrium


  !> Selects the submodule based on the specified equilibrium
  !! in the parfile. Works on a case-select basis.
  !! @warning   Throws an error if the equilibrium type is not recognised.
  subroutine set_equilibrium_pointer(settings)
    type(settings_t), intent(in) :: settings

    select case(settings%equilibrium%get_equilibrium_type())
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
    case("couette_flow")
      set_equilibrium_values => couette_flow_eq
    case("taylor_couette")
      set_equilibrium_values => taylor_couette_eq
    case("harris_sheet")
      set_equilibrium_values => harris_sheet_eq
    case("tc_pinch")
      set_equilibrium_values => tc_pinch_eq
    case("user_defined")
      set_equilibrium_values => user_defined_eq
    case default
      call logger%error( &
        "equilibrium not recognised: " &
        // trim(settings%equilibrium%get_equilibrium_type()) &
      )
    end select
  end subroutine set_equilibrium_pointer


  !> Cleaning routine, deallocates the equilibrium types.
  subroutine equilibrium_clean()
    if (allocated(grav_field%grav)) call deallocate_type(grav_field)
    if (allocated(eta_field%eta)) call deallocate_type(eta_field)
    if (allocated(rc_field%L0)) call deallocate_type(rc_field)
    if (allocated(kappa_field%kappa_para)) call deallocate_type(kappa_field)
    if (allocated(hall_field%hallfactor)) call deallocate_type(hall_field)
  end subroutine equilibrium_clean

end module mod_equilibrium
