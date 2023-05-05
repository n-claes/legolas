module mod_units
  use mod_global_variables, only: dp
  implicit none

  private

  type, public :: units_t
    logical, private :: units_set
    logical, private :: cgs
    logical, private :: based_on_density
    logical, private :: based_on_temperature
    logical, private :: based_on_numberdensity

    real(dp), private :: unit_length
    real(dp), private :: unit_time
    real(dp), private :: unit_density
    real(dp), private :: unit_velocity
    real(dp), private :: unit_temperature
    real(dp), private :: unit_pressure
    real(dp), private :: unit_magneticfield
    real(dp), private :: unit_numberdensity
    real(dp), private :: unit_mass
    real(dp), private :: mean_molecular_weight
    real(dp), private :: unit_resistivity
    real(dp), private :: unit_lambdaT
    real(dp), private :: unit_conduction

  contains

    procedure, public :: in_cgs
    procedure, public :: are_set
    procedure, public :: set_units_from_density
    procedure, public :: set_units_from_temperature
    procedure, public :: set_units_from_numberdensity
    procedure, public :: set_mean_molecular_weight
    procedure, public :: get_unit_length
    procedure, public :: get_unit_time
    procedure, public :: get_unit_density
    procedure, public :: get_unit_velocity
    procedure, public :: get_unit_temperature
    procedure, public :: get_unit_pressure
    procedure, public :: get_unit_magneticfield
    procedure, public :: get_unit_numberdensity
    procedure, public :: get_unit_mass
    procedure, public :: get_mean_molecular_weight
    procedure, public :: get_unit_resistivity
    procedure, public :: get_unit_lambdaT
    procedure, public :: get_unit_conduction
    procedure, public :: get_unit_gravity

    procedure, private :: update_dependent_units
    procedure, private :: set_based_on_to_false

  end type units_t

  public :: new_unit_system

contains

  pure function new_unit_system() result(units)
    type(units_t) :: units

    units%units_set = .false.
    units%cgs = .true.
    units%based_on_density = .false.
    units%based_on_temperature = .false.
    units%based_on_numberdensity = .false.
    units%mean_molecular_weight = 0.5_dp

    call units%set_units_from_temperature( &
      unit_length=1.0e9_dp, &
      unit_magneticfield=10.0_dp, &
      unit_temperature=1.0e6_dp &
    )
  end function new_unit_system


  pure logical function in_cgs(this)
    class(units_t), intent(in) :: this
    in_cgs = this%cgs
  end function in_cgs


  pure logical function are_set(this)
    class(units_t), intent(in) :: this
    are_set = this%units_set
  end function are_set


  pure subroutine set_based_on_to_false(this)
    class(units_t), intent(inout) :: this
    this%based_on_density = .false.
    this%based_on_temperature = .false.
    this%based_on_numberdensity = .false.
  end subroutine set_based_on_to_false


  pure subroutine set_units_from_density( &
    this, unit_length, unit_magneticfield, unit_density, mean_molecular_weight &
  )
    class(units_t), intent(inout) :: this
    real(dp), intent(in) :: unit_length
    real(dp), intent(in) :: unit_magneticfield
    real(dp), intent(in) :: unit_density
    real(dp), intent(in), optional :: mean_molecular_weight

    call this%set_based_on_to_false()
    this%unit_length = unit_length
    this%unit_magneticfield = unit_magneticfield
    this%unit_density = unit_density
    this%based_on_density = .true.
    if (present(mean_molecular_weight)) then
      this%mean_molecular_weight = mean_molecular_weight
    end if
    call this%update_dependent_units()
  end subroutine set_units_from_density


  pure subroutine set_units_from_temperature( &
    this, unit_length, unit_magneticfield, unit_temperature, mean_molecular_weight &
  )
    class(units_t), intent(inout) :: this
    real(dp), intent(in) :: unit_length
    real(dp), intent(in) :: unit_magneticfield
    real(dp), intent(in) :: unit_temperature
    real(dp), intent(in), optional :: mean_molecular_weight

    call this%set_based_on_to_false()
    this%unit_length = unit_length
    this%unit_magneticfield = unit_magneticfield
    this%unit_temperature = unit_temperature
    this%based_on_temperature = .true.
    if (present(mean_molecular_weight)) then
      this%mean_molecular_weight = mean_molecular_weight
    end if
    call this%update_dependent_units()
  end subroutine set_units_from_temperature


  pure subroutine set_units_from_numberdensity( &
    this, unit_length, unit_temperature, unit_numberdensity, mean_molecular_weight &
  )
    class(units_t), intent(inout) :: this
    real(dp), intent(in) :: unit_length
    real(dp), intent(in) :: unit_temperature
    real(dp), intent(in) :: unit_numberdensity
    real(dp), intent(in), optional :: mean_molecular_weight

    call this%set_based_on_to_false()
    this%unit_length = unit_length
    this%unit_temperature = unit_temperature
    this%unit_numberdensity = unit_numberdensity
    this%based_on_numberdensity = .true.
    if (present(mean_molecular_weight)) then
      this%mean_molecular_weight = mean_molecular_weight
    end if
    call this%update_dependent_units()
  end subroutine set_units_from_numberdensity


  pure subroutine update_dependent_units(this)
    use mod_physical_constants, only: mu0_cgs, mp_cgs, kB_cgs

    class(units_t), intent(inout) :: this

    if (this%based_on_numberdensity) then
      this%unit_density = mp_cgs * this%unit_numberdensity
      this%unit_pressure = ( &
        this%mean_molecular_weight &
        * this%unit_numberdensity &
        * kB_cgs &
        * this%unit_temperature &
      )
      this%unit_velocity = sqrt(this%unit_pressure / this%unit_density)
      this%unit_magneticfield = sqrt(mu0_cgs * this%unit_pressure)
    else if (this%based_on_density) then
      this%unit_pressure = this%unit_magneticfield**2 / mu0_cgs
      this%unit_temperature = ( &
        this%mean_molecular_weight &
        * this%unit_pressure &
        * mp_cgs &
        / (kB_cgs * this%unit_density) &
      )
      this%unit_numberdensity = this%unit_density / mp_cgs
      this%unit_velocity = this%unit_magneticfield / sqrt(mu0_cgs * this%unit_density)
    else if (this%based_on_temperature) then
      this%unit_pressure = this%unit_magneticfield**2 / mu0_cgs
      this%unit_density = ( &
        this%mean_molecular_weight &
        * this%unit_pressure &
        * mp_cgs &
        / (kB_cgs * this%unit_temperature) &
      )
      this%unit_numberdensity = this%unit_density / mp_cgs
      this%unit_velocity = this%unit_magneticfield / sqrt(mu0_cgs * this%unit_density)
    end if
    this%unit_mass = this%unit_density * this%unit_length**3
    this%unit_time = this%unit_length / this%unit_velocity
    this%unit_resistivity = this%unit_length**2 / this%unit_time
    this%unit_lambdaT = this%unit_pressure / ( &
      this%unit_time * this%unit_numberdensity**2 &
    )
    this%unit_conduction = ( &
      this%unit_density &
      * this%unit_length &
      * this%unit_velocity**3 &
      / this%unit_temperature &
    )

    this%units_set = .true.
  end subroutine update_dependent_units


  pure subroutine set_mean_molecular_weight(this, mean_molecular_weight)
    class(units_t), intent(inout) :: this
    real(dp), intent(in) :: mean_molecular_weight
    this%mean_molecular_weight = mean_molecular_weight
    if (this%units_set) call this%update_dependent_units()
  end subroutine set_mean_molecular_weight


  pure real(dp) function get_unit_length(this)
    class(units_t), intent(in) :: this
    get_unit_length = this%unit_length
  end function get_unit_length


  pure real(dp) function get_unit_time(this)
    class(units_t), intent(in) :: this
    get_unit_time = this%unit_time
  end function get_unit_time


  pure real(dp) function get_unit_density(this)
    class(units_t), intent(in) :: this
    get_unit_density = this%unit_density
  end function get_unit_density


  pure real(dp) function get_unit_velocity(this)
    class(units_t), intent(in) :: this
    get_unit_velocity = this%unit_velocity
  end function get_unit_velocity


  pure real(dp) function get_unit_temperature(this)
    class(units_t), intent(in) :: this
    get_unit_temperature = this%unit_temperature
  end function get_unit_temperature


  pure real(dp) function get_unit_pressure(this)
    class(units_t), intent(in) :: this
    get_unit_pressure = this%unit_pressure
  end function get_unit_pressure


  pure real(dp) function get_unit_magneticfield(this)
    class(units_t), intent(in) :: this
    get_unit_magneticfield = this%unit_magneticfield
  end function get_unit_magneticfield


  pure real(dp) function get_unit_numberdensity(this)
    class(units_t), intent(in) :: this
    get_unit_numberdensity = this%unit_numberdensity
  end function get_unit_numberdensity


  pure real(dp) function get_unit_mass(this)
    class(units_t), intent(in) :: this
    get_unit_mass = this%unit_mass
  end function get_unit_mass


  pure real(dp) function get_mean_molecular_weight(this)
    class(units_t), intent(in) :: this
    get_mean_molecular_weight = this%mean_molecular_weight
  end function get_mean_molecular_weight


  pure real(dp) function get_unit_resistivity(this)
    class(units_t), intent(in) :: this
    get_unit_resistivity = this%unit_resistivity
  end function get_unit_resistivity


  pure real(dp) function get_unit_lambdaT(this)
    class(units_t), intent(in) :: this
    get_unit_lambdaT = this%unit_lambdaT
  end function get_unit_lambdaT


  pure real(dp) function get_unit_conduction(this)
    class(units_t), intent(in) :: this
    get_unit_conduction = this%unit_conduction
  end function get_unit_conduction


  pure real(dp) function get_unit_gravity(this)
    class(units_t), intent(in) :: this
    get_unit_gravity = this%unit_length / this%unit_time**2
  end function get_unit_gravity
end module mod_units
