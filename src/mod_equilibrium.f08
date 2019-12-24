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
  use mod_global_variables, only: dp, gauss_gridpts, flow, radiative_cooling, thermal_conduction, resistivity, &
                                  external_gravity, geometry, k2, k3, x_start, x_end
  use mod_physical_constants, only: dpi
  use mod_grid, only: grid_gauss, initialise_grid
  implicit none

  private

  !! Default equilibrium variables
  !> Equilibrium density
  real(dp), allocatable         :: rho0_eq(:)
  !> Equilibrium temperature
  real(dp), allocatable         :: T0_eq(:)
  !> Equilibrium magnetic field in the y or theta-direction
  real(dp), allocatable         :: B02_eq(:)
  !> Equilibrium magnetic field in the z direction
  real(dp), allocatable         :: B03_eq(:)
  !> Equilibrium magnetic field
  real(dp), allocatable         :: B0_eq(:)

  !! Flow equilibrium variables
  !> Equilibrium velocity in the y or theta-direction
  real(dp), allocatable         :: v02_eq(:)
  !> Equilibrium velocity in the z direction
  real(dp), allocatable         :: v03_eq(:)

  !! Gravity equilibrium variables
  !> Equilibrium gravity
  real(dp), allocatable         :: grav_eq(:)

  !! Thermal conduction equilibrium variables
  !> Equilibrium parallel conduction
  real(dp), allocatable         :: tc_para_eq(:)
  !> Equilibrium perpendicular conduction
  real(dp), allocatable         :: tc_perp_eq(:)

  !! Radiative cooling equilibrium variables
  !> Equilibrium heat-loss function
  real(dp), allocatable         :: heat_loss_eq(:)

  !! Resistivity equilibrium variables
  !> Equilibrium resistivity
  real(dp), allocatable         :: eta_eq(:)

  public :: rho0_eq, T0_eq, B02_eq, B03_eq, B0_eq, v02_eq, v03_eq
  public :: grav_eq
  public :: tc_para_eq, tc_perp_eq, heat_loss_eq, eta_eq

  public :: initialise_equilibrium
  public :: set_equilibrium
  public :: equilibrium_clean


contains

  !> Initialises the equilibrium by allocating all equilibrium arrays and
  !! setting them to zero.
  subroutine initialise_equilibrium()
    use mod_radiative_cooling, only: initialise_radiative_cooling

    allocate(rho0_eq(gauss_gridpts))
    allocate(T0_eq(gauss_gridpts))
    allocate(B02_eq(gauss_gridpts))
    allocate(B03_eq(gauss_gridpts))
    allocate(B0_eq(gauss_gridpts))

    allocate(v02_eq(gauss_gridpts))
    allocate(v03_eq(gauss_gridpts))

    allocate(grav_eq(gauss_gridpts))

    allocate(tc_para_eq(gauss_gridpts))
    allocate(tc_perp_eq(gauss_gridpts))

    allocate(heat_loss_eq(gauss_gridpts))

    allocate(eta_eq(gauss_gridpts))

    rho0_eq = 0.0d0
    T0_eq   = 0.0d0
    B02_eq  = 0.0d0
    B03_eq  = 0.0d0
    B0_eq   = 0.0d0

    v02_eq  = 0.0d0
    v03_eq  = 0.0d0

    grav_eq = 0.0d0

    tc_para_eq = 0.0d0
    tc_perp_eq = 0.0d0

    heat_loss_eq = 0.0d0

    eta_eq = 0.0d0

    ! radiative cooling if desired
    if (radiative_cooling) then
      call initialise_radiative_cooling()
    end if

  end subroutine initialise_equilibrium

  !> Determines which pre-coded equilibrium configuration has to be loaded.
  subroutine set_equilibrium()
    use mod_global_variables, only: equilibrium_type
    use mod_check_values, only: check_negative_array
    use mod_thermal_conduction, only: get_kappa_para, get_kappa_perp
    use mod_resistivity, only: get_eta
    use mod_equilibrium_derivatives, only: set_cooling_derivatives, set_resistivity_derivatives, &
                                           set_conduction_derivatives

    select case(equilibrium_type)
      case("Adiabatic homogeneous")
        call adiabatic_homo_eq()

      case("Resistive homogeneous")
        call resistive_homo_eq()

      case("Gravitational homogeneous")
        call gravity_homo_eq()

      case("Resistive tearing modes")
        call resistive_tearing_modes_eq()

      case("Resistive tearing modes with flow")
        call resistive_tearing_modes_flow_eq()

      case("Flow driven instabilities")
        call flow_driven_instabilities_eq()

      case("Suydam cluster modes")
        call suydam_cluster_eq()

      case("Kelvin-Helmholtz")
        call KH_instability_eq()

      case("Rotating plasma cylinder")
        call rotating_plasma_cyl_eq()

      case default
        write(*, *) "Equilibrium not recognised."
        write(*, *) "Currently set on: ", equilibrium_type
        stop
    end select

    !! Set equilibrium physics
    if (radiative_cooling) then
      call set_cooling_derivatives(T0_eq, rho0_eq)
      ! this is L_0, should balance out in thermal equilibrium.
      heat_loss_eq = 0.0d0
    end if
    if (thermal_conduction) then
      call get_kappa_para(T0_eq, tc_para_eq)
      call get_kappa_perp(T0_eq, rho0_eq, B0_eq, tc_perp_eq)
      call set_conduction_derivatives(T0_eq, rho0_eq, B0_eq)
    end if
    if (resistivity) then
      call get_eta(T0_eq, eta_eq)
      call set_resistivity_derivatives(T0_eq)
    end if

    ! Check for negative pressure, temperature, etc.
    call check_negative_array(rho0_eq, "density")
    call check_negative_array(T0_eq, "temperature")
  end subroutine set_equilibrium


  !> Initialises equilibrium for an adiabatic homogeneous medium in Cartesian
  !! geometry.
  subroutine adiabatic_homo_eq()
    geometry = 'Cartesian'
    call initialise_grid()

    flow = .false.
    radiative_cooling = .false.
    thermal_conduction = .false.
    resistivity = .false.
    external_gravity = .false.

    k2 = dpi
    k3 = dpi

    !! Parameters
    rho0_eq = 1.0d0
    T0_eq   = 1.0d0
    B02_eq  = 1.0d0
    B03_eq  = 1.0d0
    B0_eq   = sqrt(B02_eq**2 + B03_eq**2)
  end subroutine adiabatic_homo_eq


  !> Initialises equilibrium for a homogeneous medium in Cartesian geometry
  !! with a constant resistivity value. From Advanced MHD, page 156.
  subroutine resistive_homo_eq()
    use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value

    real(dp)  :: beta

    geometry = 'Cartesian'
    call initialise_grid()

    flow = .false.
    radiative_cooling = .false.
    thermal_conduction = .false.
    resistivity = .true.
    use_fixed_resistivity = .true.
    fixed_eta_value = 0.001d0
    external_gravity = .false.

    k2 = 0.0d0
    k3 = 1.0d0
    beta = 0.25d0

    !! Parameters
    rho0_eq = 1.0d0
    B02_eq  = 0.0d0
    B03_eq  = 1.0d0
    B0_eq   = sqrt(B02_eq**2 + B03_eq**2)
    T0_eq   = beta * B0_eq**2 / (2) ! n=1, kB=1, mu0=1
  end subroutine resistive_homo_eq


  !> Initialises equilibrium for a homogeneous medium in Cartesian geometry
  !! with a constant gravity term included.
  subroutine gravity_homo_eq()
    geometry = 'Cartesian'
    call initialise_grid()

    flow = .false.
    radiative_cooling = .false.
    thermal_conduction = .false.
    resistivity = .false.
    external_gravity = .true.

    k2 = dpi
    k3 = dpi

    !! Parameters
    rho0_eq = 1.0d0
    B02_eq  = 1.0d0
    B03_eq  = 1.0d0
    B0_eq   = sqrt(B02_eq**2 + B03_eq**2)
    T0_eq   = 1.0d0
    grav_eq = 0.5d0
  end subroutine gravity_homo_eq


  !> Initialises equilibrium for the resistive tearing modes in Cartesian
  !! geometry, without flow. From Advanced MHD, page 159.
  subroutine resistive_tearing_modes_eq()
    use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value
    use mod_equilibrium_derivatives, only: d_B02_dr, d_B03_dr, dd_B02_dr, dd_B03_dr

    real(dp)              :: alpha, beta, x
    integer               :: i

    geometry = 'Cartesian'
    ! Override values from par file
    x_start = -0.5d0
    x_end   = 0.5d0
    call initialise_grid()

    flow = .false.
    radiative_cooling = .false.
    thermal_conduction = .false.
    resistivity = .true.
    use_fixed_resistivity = .true.
    fixed_eta_value = 0.0001d0
    external_gravity = .false.

    k2 = 0.49d0
    k3 = 0.0d0

    ! Parameters
    alpha = 4.73884d0
    beta  = 0.15d0

    do i = 1, gauss_gridpts
      x = grid_gauss(i)

      ! Equilibrium
      rho0_eq(i) = 1.0d0
      B02_eq(i)  = sin(alpha * x)
      B03_eq(i)  = cos(alpha * x)
      B0_eq(i)   = sqrt(B02_eq(i)**2 + B03_eq(i)**2)
      T0_eq(i)   = beta * B0_eq(i)**2 / (2) ! n=1, kB=1, mu0=1

      ! Derivatives
      d_B02_dr(i)   = alpha * cos(alpha * x)
      d_B03_dr(i)   = -alpha * sin(alpha * x)
      ! No d_T0_dr needed, as B0**2 is independent of r

      dd_B02_dr(i)  = -alpha**2 * sin(alpha * x)
      dd_B03_dr(i)  = -alpha**2 * cos(alpha * x)
    end do
  end subroutine resistive_tearing_modes_eq


  !> Initialises equilibrium for the resistive tearing modes in Cartesian
  !! geometry, with flow. From Advanced MHD, page 161.
  subroutine resistive_tearing_modes_flow_eq()
    use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value
    use mod_equilibrium_derivatives, only: d_B02_dr, d_B03_dr, dd_B02_dr, dd_B03_dr, d_v02_dr

    real(dp)    :: alpha, beta, x
    integer     :: i

    geometry = 'Cartesian'
    ! Override values from par file
    x_start = -0.5d0
    x_end   = 0.5d0
    call initialise_grid()

    flow = .true.
    radiative_cooling = .false.
    thermal_conduction = .false.
    resistivity = .true.
    use_fixed_resistivity = .true.
    fixed_eta_value = 0.0001d0
    external_gravity = .false.

    k2 = 1.5d0
    k3 = 0.0d0

    ! Parameters
    alpha = 4.73884d0
    beta  = 0.15d0

    do i = 1, gauss_gridpts
      x = grid_gauss(i)

      ! Equilibrium
      rho0_eq(i) = 1.0d0
      B02_eq(i)  = sin(alpha * x)
      B03_eq(i)  = cos(alpha * x)
      B0_eq(i)   = sqrt(B02_eq(i)**2 + B03_eq(i)**2)
      T0_eq(i)   = beta * B0_eq(i)**2 / (2) ! n=1, kB=1, mu0=1
      v02_eq(i)  = 0.15 * x

      ! Derivatives
      d_B02_dr(i)   = alpha * cos(alpha * x)
      d_B03_dr(i)   = -alpha * sin(alpha * x)
      ! No d_T0_dr needed, as B0**2 is independent of r

      dd_B02_dr(i)  = -alpha**2 * sin(alpha * x)
      dd_B03_dr(i)  = -alpha**2 * cos(alpha * x)

      d_v02_dr(i)  = 0.15
    end do
  end subroutine resistive_tearing_modes_flow_eq


  subroutine flow_driven_instabilities_eq()
    use mod_equilibrium_derivatives, only: d_rho0_dr, d_v02_dr, d_v03_dr, d_B02_dr, d_B03_dr, d_T0_dr

    real(dp)    :: rho0, delta, theta, v0, v1, v2, tau, phi0, alpha, B0, x, k0, g, p0
    real(dp)    :: v_x(gauss_gridpts), phi_x(gauss_gridpts), p_x(gauss_gridpts)
    integer     :: i

    geometry = 'Cartesian'
    call initialise_grid()

    g = 15.0d0

    flow = .true.
    radiative_cooling = .false.
    thermal_conduction = .false.
    resistivity = .false.
    external_gravity = .true.

    ! Parameters
    k0 = 1.0d0
    delta = -5.0d0
    phi0 = -0.35d0 * dpi
    alpha = 0.0d0
    theta = 0.35d0 * dpi
    v0 = 0.2d0
    v1 = 0.6d0
    v2 = 0.0d0
    tau = 0.0d0

    k2 = 0.0d0
    k3 = k0

    rho0 = 1.0d0
    p0   = 1.0d0
    B0   = 1.0d0
    grav_eq = g

    do i = 1, gauss_gridpts
      x = grid_gauss(i)
      v_x(i)   = v0 + v1*(x - 0.5d0) + v2*sin(tau * (x - 0.5d0))
      phi_x(i) = phi0 + alpha*(x - 0.5d0)
      p_x(i)   = p0 - (x - 0.5d0 * delta * x**2)*g

      rho0_eq(i) = rho0 * (1 - delta*x)
      v02_eq(i)  = sin(theta) * v_x(i)
      v03_eq(i)  = cos(theta) * v_x(i)
      B02_eq(i)  = B0 * sin(phi_x(i))
      B03_eq(i)  = B0 * cos(phi_x(i))
      B0_eq(i)   = sqrt(B02_eq(i)**2 + B03_eq(i)**2)
      T0_eq(i)   = p_x(i) / rho0_eq(i)

      d_rho0_dr(i) = -rho0 * delta
      d_v02_dr(i)  = sin(theta) * (v1 + v2*cos(tau * (x - 0.5d0)) * tau)
      d_v03_dr(i)  = cos(theta) * (v1 + v2*cos(tau * (x - 0.5d0)) * tau)
      d_B02_dr(i)  = B0 * cos(phi_x(i)) * alpha
      d_B03_dr(i)  = -B0 * sin(phi_x(i)) * alpha
      d_T0_dr(i)   = (-rho0*(1.0d0 - delta * x)**2 * g + rho0*delta*p_x(i)) &
                      / (rho0 * delta)**2
    end do
  end subroutine flow_driven_instabilities_eq


  !> Initialises equilibrium for Suydam cluster modes in cylindrical geometry.
  !! Obtained from Bondeson et al., Phys. Fluids 30 (1987)
  subroutine suydam_cluster_eq()
    use mod_equilibrium_derivatives, only: d_T0_dr, d_v03_dr

    real(dp)      :: v_z0, p0, p1, alpha, r
    real(dp)      :: J0, J1, DJ0, DJ1
    real(dp)      :: P0_eq(gauss_gridpts)
    integer       :: i

    geometry = 'cylindrical'
    call initialise_grid()

    flow = .true.
    radiative_cooling = .false.
    thermal_conduction = .false.
    resistivity = .false.
    external_gravity = .false.

    !! Parameters
    v_z0  = 0.14d0
    p0    = 0.05d0
    p1    = 0.1d0
    alpha = 2.0d0

    k2 = 1.0d0
    k3 = -1.2d0

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      J0 = bessel_jn(0, alpha * r)
      J1 = bessel_jn(1, alpha * r)
      DJ0 = -alpha * J1
      DJ1 = alpha * (0.5d0 * J0 - 0.5d0 * bessel_jn(2, alpha * r))

      ! Equilibrium
      rho0_eq(i) = 1.0d0
      v02_eq(i) = 0.0d0
      v03_eq(i) = v_z0 * (1.0d0 - r**2)
      B02_eq(i) = J1
      B03_eq(i) = sqrt(1.0d0 - p1) * J0
      B0_eq(i)  = sqrt(B02_eq(i)**2 + B03_eq(i)**2)
      P0_eq(i)  = p0 + 0.5d0 * p1 * J0**2
      T0_eq(i)  = P0_eq(i) / rho0_eq(i)

      ! Derivatives
      d_T0_dr(i)    = p1 * J0 * DJ0
      d_v03_dr(i)   = -2.0d0 * v_z0 * r
    end do
  end subroutine suydam_cluster_eq


  !> Initialises equilibrium for the Kelvin-Helmholtz instability in
  !! Cartesian geometry.
  !! From Miura et al., J. Geophys. Res. 87 (1982)
  subroutine KH_instability_eq()
    use mod_equilibrium_derivatives, only: d_v02_dr, d_v03_dr

    real(dp)    :: a, x, p0, v0y, v0z
    integer     :: i

    geometry = 'Cartesian'
    ! Override values from par file: if on interval [0,1], change velocity profile x --> x - 0.5
    x_start = -0.5d0
    x_end   = 0.5d0
    call initialise_grid()

    flow = .true.
    radiative_cooling = .false.
    thermal_conduction = .false.
    resistivity = .false.
    external_gravity = .true.

    !! Parameters
    a   = 0.05d0
    P0  = 3.6d0
    v0y = 1.67d0
    v0z = 0.0d0   ! so B is perpendicular to v

    !! Filling equilibria
    rho0_eq = 1.0d0
    B02_eq  = 0.0d0
    B03_eq  = 1.0d0
    B0_eq   = sqrt(B02_eq**2 + B03_eq**2)
    grav_eq = 10.0d0

    k2 = 10.0d0
    k3 = 0.0d0

    do i = 1, gauss_gridpts
      x = grid_gauss(i)
      v02_eq(i)   = -v0y * tanh(x / a)
      v03_eq(i)   = -v0z * tanh(x / a)

      T0_eq(i)    = P0 / rho0_eq(i)

      ! Derivatives
      d_v02_dr(i) = -v0y / (a * cosh(x / a)**2)
      d_v03_dr(i) = -v0z / (a * cosh(x / a)**2)
    end do
  end subroutine KH_instability_eq


  subroutine rotating_plasma_cyl_eq()
    use mod_equilibrium_derivatives, only: d_B02_dr, d_B03_dr, d_v02_dr, d_v03_dr, d_T0_dr

    real(dp)    :: a21, a22, a3, b21, b22, b3, p0
    real(dp)    :: r
    integer     :: i

    geometry = 'cylindrical'
    call initialise_grid() ! Initialise grid

    flow = .true.
    radiative_cooling = .false.
    thermal_conduction = .false.
    resistivity = .false.
    external_gravity = .false.

    !! Parameters
    a21 = 8.0d0
    a22 = 0.0d0
    a3  = 0.0d0
    b21 = 1.0d0
    b22 = 0.0d0
    b3  = 0.0d0
    p0  = 0.1d0

    k2  = 1.0d0
    k3  = 0.0d0

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      !! Equilibrium
      rho0_eq(i) = 1.0d0
      v02_eq(i)  = a21*r + a22*r**2
      v03_eq(i)  = a3
      B02_eq(i)  = b21*r + b22*r**2
      B03_eq(i)  = b3
      B0_eq(i)   = sqrt(B02_eq(i)**2 + B03_eq(i)**2)
      T0_eq(i)   = (1.0d0 / rho0_eq(i)) * (p0 &
                      + 0.5d0 * (a21**2 - 2*b21**2)*r**2 &
                      + (2.0d0/3.0d0)*(a21*a22 - b21*b22)*r**3 &
                      + (1.0d0/4.0d0)*(a22**2 - b22**2)*r**4)

      !! Derivatives
      d_B02_dr(i) = b21 + 2.0d0*b22*r
      d_B03_dr(i) = 0.0d0
      d_v02_dr(i) = a21 + 2.0d0*a22*r
      d_v03_dr(i) = 0.0d0
      d_T0_dr(i)  = (1.0d0 / rho0_eq(i)) * ( (a21**2 - 2.0d0*b21**2)*r &
                     + 2.0d0*(a21*a22 - b21*b22)*r**2 &
                     + (a22**2 - b22**2)*r**3 )
    end do
  end subroutine rotating_plasma_cyl_eq


  !> Cleaning routine, deallocates all arrays in this module.
  subroutine equilibrium_clean()
    deallocate(rho0_eq)
    deallocate(T0_eq)
    deallocate(B02_eq)
    deallocate(B03_eq)
    deallocate(B0_eq)
    deallocate(v02_eq)
    deallocate(v03_eq)
    deallocate(grav_eq)
    deallocate(tc_para_eq)
    deallocate(tc_perp_eq)
    deallocate(heat_loss_eq)
    deallocate(eta_eq)
  end subroutine equilibrium_clean


end module mod_equilibrium
