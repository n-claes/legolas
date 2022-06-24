! =============================================================================
!> Module defining all equilibrium types and governing initialisation.
!! This module contains all equilibrium types and contains the function
!! allocating and deallocating them.
module mod_equilibrium_fields
  use mod_types
  implicit none

  private

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

end module mod_equilibrium_fields
