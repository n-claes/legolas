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
      procedure(profile_fcn), pointer, nopass :: v   => null()
      procedure(profile_fcn), pointer, nopass :: dv  => null()
    end type ic_velocity_t

    type, public :: ic_temperature_t
      procedure(profile_fcn), pointer, nopass :: T   => null()
      procedure(profile_fcn), pointer, nopass :: dT  => null()
    end type ic_temperature_t
  
    type, public :: initial_conditions_t
      type(ic_density_t)     :: density
      type(ic_velocity_t)    :: velocity_1
      type(ic_velocity_t)    :: velocity_2
      type(ic_velocity_t)    :: velocity_3
      type(ic_temperature_t) :: temperature
    contains
      procedure :: set_ic_density_funcs
      procedure :: set_ic_velocity_1_funcs
      procedure :: set_ic_velocity_2_funcs
      procedure :: set_ic_velocity_3_funcs
      procedure :: set_ic_temperature_funcs
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
  
      ic%velocity_1%v  => zero_fcn
      ic%velocity_1%dv => zero_fcn

      ic%velocity_2%v  => zero_fcn
      ic%velocity_2%dv => zero_fcn

      ic%velocity_3%v  => zero_fcn
      ic%velocity_3%dv => zero_fcn

      ic%temperature%T  => zero_fcn
      ic%temperature%dT => zero_fcn

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
      procedure(profile_fcn) :: dv01_func
  
      call logger%debug("Setting ICs for component v1.")
      self%velocity_1%v => v01_func
      self%velocity_1%dv => dv01_func
    end subroutine set_ic_velocity_1_funcs

    subroutine set_ic_velocity_2_funcs(self, v02_func, dv02_func)
      class(initial_conditions_t), intent(inout) :: self
      procedure(profile_fcn) :: v02_func
      procedure(profile_fcn) :: dv02_func

      call logger%debug("Setting ICs for component v2.")
      self%velocity_2%v => v02_func
      self%velocity_2%dv => dv02_func
    end subroutine set_ic_velocity_2_funcs

    subroutine set_ic_velocity_3_funcs(self, v03_func, dv03_func)
      class(initial_conditions_t), intent(inout) :: self
      procedure(profile_fcn) :: v03_func
      procedure(profile_fcn) :: dv03_func

      call logger%debug("Setting ICs for component v3.")
      self%velocity_3%v => v03_func
      self%velocity_3%dv => dv03_func
    end subroutine set_ic_velocity_3_funcs

    subroutine set_ic_temperature_funcs(self, T_func, dT_func)
      class(initial_conditions_t), intent(inout) :: self
      procedure(profile_fcn) :: T_func
      procedure(profile_fcn), optional :: dT_func

      call logger%debug("Setting ICs for component T.")
      self%temperature%T => T_func
      if (present(dT_func)) self%temperature%dT => dT_func
    end subroutine set_ic_temperature_funcs
  
  end module mod_iv_initial_conditions
  