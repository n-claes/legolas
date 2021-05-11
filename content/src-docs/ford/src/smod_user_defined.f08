!
! SUBMODULE: smod_user_defined
!
! DESCRIPTION:
! Submodule for user defined equilibria or tests
submodule (mod_equilibrium) smod_user_defined
  implicit none

contains

  module subroutine user_defined_eq()
    !To uncomment in case of a fixed resistivity
    !use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value
    ! To uncomment in case of a fixed thermal conduction
    !use mod_global_variables, only: use_fixed_tc, fixed_tc_para_value, fixed_tc_perp_value

    !! Declare additional parameters
    real(dp)    :: x, alpha, beta, v0
    integer     :: i

    ! 'Cartesian' or 'cylindrical' (caps sensitive)
    geometry = 'Cartesian'
    x_start = 0.0d0
    x_end   = 1.0d0
    call initialise_grid()

    ! Wave numbers:
    ! Cartesian:   k2 = k_y ; k3 = k_z
    ! cylindrical: k2 = m   ; k3 = k
    k2 = 1.5d0
    k3 = 0.0d0

    ! Uncomment in case of flow
    !flow = .true.

    ! Uncomment in case of resistivity
    !resistivity = .true.
    ! Uncomment when using fixed resistivity and define value
    !use_fixed_resistivity = .true.
    !fixed_eta_value = 0.0001d0

    ! Uncomment in case of thermal conduction
    !thermal_conduction = .true.
    ! Uncomment when using fixed thermal conduction and define values
    !use_fixed_tc = .true.
    !fixed_tc_para_value = 0.001d0
    !fixed_tc_perp_value = 0.0d0

    ! Uncomment in case of radiative cooling
    !radiative_cooling = .true.

    ! Uncomment in case of external gravity
    !external_gravity = .true.

    ! Define parameters
    alpha = 4.73884d0
    beta  = 0.15d0
    v0    = 0.15d0

    do i = 1, gauss_gridpts
      x = grid_gauss(i)

      ! Equilibrium
      rho_field % rho0(i) = 1.0d0
      B_field % B02(i)    = sin(alpha * x)
      B_field % B03(i)    = cos(alpha * x)
      B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      T_field % T0(i)     = beta * (B_field % B0(i))**2 / (2.0d0)

      ! Derivatives
      rho_field % d_rho0_dr(i) = 0.0d0
      B_field % d_B02_dr(i)    = alpha * cos(alpha * x)
      B_field % d_B03_dr(i)    = -alpha * sin(alpha * x)
      T_field % d_T0_dr(i)     = 0.0d0

      ! Define in case of flow
      !v_field % v02(i)      = v0 * sin(x)
      !v_field % v03(i)      = v0 * cos(x)
      !v_field % d_v02_dr(i) = v0

      ! Define in case of resistivity
      !eta_field % dd_B02_dr(i) = -alpha**2 * sin(alpha * x)
      !eta_field % dd_B03_dr(i) = -alpha**2 * cos(alpha * x)

      ! Define in case of gravity
      !grav_field % grav(i) = 0.5d0
    end do

  end subroutine user_defined_eq

end submodule smod_user_defined
