module mod_iv_initial_conditions
    use mod_global_variables, only: dp
    use mod_logging,          only: logger
    use mod_iv_globals,       only: iv_fcn_ptr_t, profile_fcn, zero_fcn
    implicit none
  
    private
  
    type, public :: ic_density_t
      procedure(profile_fcn), pointer, nopass :: rho   => null()
      procedure(profile_fcn), pointer, nopass :: drho  => null()
    end type ic_density_t
  
    type, public :: ic_velocity_t
      procedure(profile_fcn), pointer, nopass :: v01   => null()
      procedure(profile_fcn), pointer, nopass :: dv01  => null()
      ! v02, etc.
    end type ic_velocity_t
  
    type, public :: initial_conditions_t
      type(ic_density_t)   :: density
      type(ic_velocity_t)  :: velocity_1
      ! also temperature, etc.
    contains
      procedure :: set_ic_density_funcs
      procedure :: set_ic_velocity_1_funcs
    end type initial_conditions_t

    public :: new_initial_conditions
  
  contains
  
    !=====================================================
    ! Constructor
    !=====================================================
    function new_initial_conditions() result(ic)
      type(initial_conditions_t) :: ic
      
      ! Default everything to zero
      ic%density%rho   => zero_fcn
      ic%density%drho  => zero_fcn
  
      ic%velocity_1%v01  => zero_fcn
      ic%velocity_1%dv01 => zero_fcn
  
      ! TODO: add the rest

    end function new_initial_conditions
  
    !=====================================================
    ! Setter routines
    !=====================================================
    subroutine set_ic_density_funcs(self, rho_func, drho_func)
      class(initial_conditions_t), intent(inout) :: self
      procedure(profile_fcn) :: rho_func
      procedure(profile_fcn), optional :: drho_func
  
      call logger%debug("Setting ICs for component rho.")
      self%density%rho => rho_func
      if (present(drho_func)) self%density%drho => drho_func
    end subroutine set_ic_density_funcs

  
    subroutine set_ic_velocity_1_funcs(self, v01_func, dv01_func)
      class(initial_conditions_t), intent(inout) :: self
      procedure(profile_fcn) :: v01_func
      procedure(profile_fcn), optional :: dv01_func
  
      call logger%debug("Setting ICs for component v1.")
      self%velocity_1%v01 => v01_func
      if (present(dv01_func)) self%velocity_1%dv01 => dv01_func
    end subroutine set_ic_velocity_1_funcs
  
    ! TODO: Add as needed
  
  end module mod_iv_initial_conditions
  