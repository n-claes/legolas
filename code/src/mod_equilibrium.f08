module mod_equilibrium
  use mod_global_variables
  implicit none

  !> Wavenumber in y-direction (Cartesian) or theta-direction (cylindrical)
  real(dp)                      :: k2
  !> Wavenumber in z-direction (Cartesian) or z-direction (cylindrical)
  real(dp)                      :: k3

  !! Default equilibrium variables
  !> equilibrium density
  real(dp), allocatable         :: rho0_eq(:)
  !> equilibrium temperature
  real(dp), allocatable         :: T0_eq(:)
  !> equilibrium magnetic fields
  real(dp), allocatable         :: B01_eq(:), B02_eq(:), B03_eq(:), B0_eq(:)

  !! Flow equilibrium variables
  !> equilibrium velocities
  real(dp), allocatable         :: v01_eq(:), v02_eq(:), v03_eq(:)

  !! Thermal conduction equilibrium variables
  !> equilibrium parallel conduction
  real(dp), allocatable         :: tc_para_eq(:)
  !> equilibrium perpendicular conduction
  real(dp), allocatable         :: tc_perp_eq(:)

  !! Radiative cooling equilibrium variables
  !> equilibrium heat-loss function
  real(dp), allocatable         :: heat_loss_eq(:)

  !! Resistivity equilibrium variables
  real(dp), allocatable         :: eta_eq(:)



contains

  subroutine initialise_equilibrium()
    allocate(rho0_eq(4*gridpts))
    allocate(T0_eq(4*gridpts))
    allocate(B01_eq(4*gridpts))
    allocate(B02_eq(4*gridpts))
    allocate(B03_eq(4*gridpts))
    allocate(B0_eq(4*gridpts))

    allocate(v01_eq(4*gridpts))
    allocate(v02_eq(4*gridpts))
    allocate(v03_eq(4*gridpts))

    allocate(tc_para_eq(4*gridpts))
    allocate(tc_perp_eq(4*gridpts))

    allocate(heat_loss_eq(4*gridpts))

    allocate(eta_eq(4*gridpts))


    rho0_eq = 0.0d0
    T0_eq   = 0.0d0
    B01_eq  = 0.0d0
    B02_eq  = 0.0d0
    B03_eq  = 0.0d0
    B0_eq   = 0.0d0

    v01_eq  = 0.0d0
    v02_eq  = 0.0d0
    v03_eq  = 0.0d0

    tc_para_eq = 0.0d0
    tc_perp_eq = 0.0d0

    heat_loss_eq = 0.0d0

    eta_eq = 0.0d0

    k2 = 1.0d0
    k3 = 1.0d0

    call set_equilibrium()

  end subroutine initialise_equilibrium

  !> Sets the equilibrium on the nodes of the Gaussian quadrature
  subroutine set_equilibrium()
    use mod_grid
    use mod_physical_constants
    use mod_radiative_cooling
    use mod_thermal_conduction

    ! Initialisations
    rho0_eq = 1.0d0
    T0_eq   = 1.0d0
    B01_eq  = 1.0d0
    B02_eq  = 1.0d0
    B03_eq  = 1.0d0
    B0_eq   = sqrt(B01_eq**2 + B02_eq**2 + B03_eq**2)

    if (flow) then
      v01_eq  = 1.0d0
      v02_eq  = 1.0d0
      v03_eq  = 1.0d0
    end if

    if (thermal_conduction) then
      call get_kappa_para(T0_eq, tc_para_eq)
      call get_kappa_perp(T0_eq, rho0_eq, B0_eq, tc_perp_eq)
    end if

    if (radiative_cooling) then
      ! this is L_0, should balance out in thermal equilibrium.
      heat_loss_eq = 0.0d0
    end if

    if (resistivity) then
      eta_eq = 1.0d0
    end if

  end subroutine set_equilibrium

  subroutine equilibrium_clean()
    deallocate(rho0_eq)
    deallocate(T0_eq)
    deallocate(B01_eq)
    deallocate(B02_eq)
    deallocate(B03_eq)
    deallocate(B0_eq)

    deallocate(v01_eq)
    deallocate(v02_eq)
    deallocate(v03_eq)

    deallocate(tc_para_eq)
    deallocate(tc_perp_eq)

    deallocate(heat_loss_eq)

    deallocate(eta_eq)

  end subroutine equilibrium_clean


end module mod_equilibrium
