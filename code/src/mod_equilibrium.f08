module mod_equilibrium
  use mod_global_variables
  implicit none

  !> Wavenumber in y-direction (Cartesian) or theta-direction (cylindrical)
  real(dp)                  :: k2
  !> Wavenumber in z-direction (Cartesian) or z-direction (cylindrical)
  real(dp)                  :: k3
  real(dp), allocatable         :: rho0_eq(:)
  real(dp), allocatable         :: v01_eq(:), v02_eq(:), v03_eq(:)
  real(dp), allocatable         :: T0_eq(:)
  real(dp), allocatable         :: B01_eq(:), B02_eq(:), B03_eq(:), B0_eq(:)

  real(dp), allocatable         :: heat_loss_eq(:)

  real(dp), allocatable         :: tc_para_eq(:)
  real(dp), allocatable         :: tc_perp_eq(:)


contains

  subroutine initialise_equilibrium()

    allocate(rho0_eq(4*gridpts))
    allocate(v01_eq(4*gridpts))
    allocate(v02_eq(4*gridpts))
    allocate(v03_eq(4*gridpts))
    allocate(T0_eq(4*gridpts))
    allocate(B01_eq(4*gridpts))
    allocate(B02_eq(4*gridpts))
    allocate(B03_eq(4*gridpts))
    allocate(B0_eq(4*gridpts))

    allocate(heat_loss_eq(4*gridpts))

    allocate(tc_para_eq(4*gridpts))
    allocate(tc_perp_eq(4*gridpts))

    rho0_eq = 0.0d0
    v01_eq  = 0.0d0
    v02_eq  = 0.0d0
    v03_eq  = 0.0d0
    T0_eq   = 0.0d0
    B01_eq  = 0.0d0
    B02_eq  = 0.0d0
    B03_eq  = 0.0d0
    B0_eq   = 0.0d0

    heat_loss_eq = 0.0d0

    tc_para_eq = 0.0d0
    tc_perp_eq = 0.0d0

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

    real(dp)              :: B0_eq(4*gridpts)
    real(dp)              :: x_lo, x_hi, dx, xi(n_gauss)
    integer               :: i, j, idx

    ! Check for origin in cylindrical coordinates
    if (geometry == "cylindrical") then
      grid_gauss(1) = 1.0d-5
    end if

    ! Evaluate grid_gauss in nodes of Gaussian quadrature.
    ! An integral of f(x) in [a, b] can be approximated by
    ! 0.5*(b-a) * SUM[i from 1 -> n] ( wi * f( 0.5*(b-a)*xi + 0.5*(b-a)) )
    ! where wi and xi are the weights and nodes at i (so 1 to 4 here).
    ! Hence we need the gridpoints equal to the evaluation points of f(x).
    do i = 1, gridpts - 1
      x_lo = grid(i)
      x_hi = grid(i + 1)
      dx   = x_hi - x_lo

      do j = 1, n_gauss
        xi(j) = 0.5 * dx * gaussian_nodes(j) + 0.5*(x_lo + x_hi)
        idx   = (i - 1)*n_gauss + j
        grid_gauss(idx) = xi(j)
      end do
    end do

    ! Temporary initialisations
    rho0_eq = 1.0d0
    v01_eq  = 1.0d0
    v02_eq  = 1.0d0
    v03_eq  = 1.0d0
    T0_eq   = 1.0d0
    B01_eq  = 1.0d0
    B02_eq  = 1.0d0
    B03_eq  = 1.0d0
    B0_eq   = sqrt(B01_eq**2 + B02_eq**2 + B03_eq**2)

    if (radiative_cooling) then
      ! this is L_0, should balance out in thermal equilibrium.
      heat_loss_eq = 0.0d0
    end if

    if (thermal_conduction) then
      B0_eq = sqrt(B01_eq**2 + B02_eq**2 + B03_eq**2)
      call get_kappa_para(T0_eq, tc_para_eq)
      call get_kappa_perp(T0_eq, rho0_eq, B0_eq, tc_perp_eq)
    end if

  end subroutine set_equilibrium

  subroutine equilibrium_clean()
    deallocate(rho0_eq)
    deallocate(v01_eq)
    deallocate(v02_eq)
    deallocate(v03_eq)
    deallocate(T0_eq)
    deallocate(B01_eq)
    deallocate(B02_eq)
    deallocate(B03_eq)
    deallocate(B0_eq)

    deallocate(heat_loss_eq)

    deallocate(tc_para_eq)
    deallocate(tc_perp_eq)
  end subroutine equilibrium_clean


end module mod_equilibrium
