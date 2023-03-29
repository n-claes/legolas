! =============================================================================
!> Module to inspect if certain conditions are fulfilled by doing
!! additional sanity checks on the equilibrium configuration.
!! For cylindrical geometries we check if \(k_2\) is an integer and if the
!! on-axis values obey regularity conditions. Equilibrium balance
!! for both the Cartesian and cylindrical cases is checked.
module mod_inspections
  use mod_global_variables, only: dp, dp_LIMIT
  use mod_types, only: gravity_type, cooling_type, conduction_type
  use mod_logging, only: logger, str, exp_fmt
  use mod_settings, only: settings_t
  use mod_background, only: background_t
  use mod_function_utils, only: from_function
  use mod_grid, only: grid_gauss, eps_grid, d_eps_grid_dr
  use mod_check_values, only: is_NaN, is_negative, is_zero
  implicit none

  private

  public :: perform_NaN_and_negative_checks
  public :: perform_sanity_checks
  public :: check_wavenumbers

contains

  !> General routine to do initial sanity checks on the various equilibrium attributes.
  !! We check the equilibrium arrays for NaN and see if all density and temperature
  !! values are positive.
  subroutine perform_NaN_and_negative_checks(background, grav_field)
    type(background_t), intent(in) :: background
    !> the type containing the gravity attributes
    type(gravity_type), intent(in)      :: grav_field

    character(50) :: name

    name = ""
    ! TODO: is there an easier way to do this?
    if (any(is_negative(from_function(background%density%rho0, grid_gauss)))) then
      call logger%error("negative density encountered!")
      return
    end if
    if (any(is_negative(from_function(background%temperature%T0, grid_gauss)))) then
      call logger%error("negative temperature encountered!")
      return
    end if
    if (any(is_NaN(from_function(background%density%rho0, grid_gauss)))) then
      name = "density"
    else if (any(is_NaN(from_function(background%temperature%T0, grid_gauss)))) then
      name = "temperature"
    else if (any(is_NaN(from_function(background%magnetic%B01, grid_gauss)))) then
      name = "B01"
    else if (any(is_NaN(from_function(background%magnetic%B02, grid_gauss)))) then
      name = "B02"
    else if (any(is_NaN(from_function(background%magnetic%B03, grid_gauss)))) then
      name = "B03"
    else if (any(is_NaN(from_function(background%velocity%v01, grid_gauss)))) then
      name = "v01"
    else if (any(is_NaN(from_function(background%velocity%v02, grid_gauss)))) then
      name = "v02"
    else if (any(is_NaN(from_function(background%velocity%v03, grid_gauss)))) then
      name = "v03"
    else if (any(is_NaN(grav_field % grav))) then
      name = "gravity"
    end if

    if (name /= "") then
      call logger%error("NaN encountered in " // adjustl(trim(name)))
      return
    end if
  end subroutine perform_NaN_and_negative_checks


  !> General routine to do sanity checks on the different equilibrium types.
  !! We check the wavenumbers and on-axis values, as well as standard
  !! and non-adiabatic equilibrium force balance.
  subroutine perform_sanity_checks( &
    settings, background, grav_field, rc_field, kappa_field &
  )
    !> the settings object
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    !> the type containing the gravity attributes
    type(gravity_type), intent(in)      :: grav_field
    !> the type containing the radiative cooling attributes
    type(cooling_type), intent(inout)   :: rc_field
    !> the type containing the thermal conduction attributes
    type(conduction_type), intent(in)   :: kappa_field

    call check_wavenumbers(geometry=settings%grid%get_geometry())
    call check_B01_cylindrical(settings, background)
    call check_on_axis_values(settings, background)
    call standard_equil_conditions(settings, background, grav_field)
    call continuity_equil_conditions(settings, background)
    call induction_equil_conditions(settings, background)
    call check_energy_balance(settings, background, rc_field, kappa_field)
  end subroutine perform_sanity_checks


  subroutine check_B01_cylindrical(settings, background)
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background

    if (.not. settings%grid%get_geometry() == "cylindrical") return
    if (.not. all(is_zero(from_function(background%magnetic%B01, grid_gauss)))) then
      call logger%error( &
        "B01 component currently not supported for cylindrical geometries!" &
      )
    end if
  end subroutine check_B01_cylindrical


  !> Sanity check on the wavenumbers.
  !! Checks if k2 is an integer in cylindrical geometry.
  !! @warning An error if thrown if the geometry is cylindrical and k2 is not
  !!          an integer.
  subroutine check_wavenumbers(geometry)
    use mod_equilibrium_params, only: k2

    character(len=*), intent(in) :: geometry

    ! in cylindrical geometry k2 should be an integer
    if (geometry == "cylindrical" .and. abs(int(k2) - k2) > dp_LIMIT) then
      call logger%error( &
        "cylindrical geometry but k2 is not an integer! Value: " // str(k2) &
      )
    end if
  end subroutine check_wavenumbers


  !> Checks if on-axis regularity conditions are satisfied in cylindrical geometry.
  !! We check if \(B_\theta, B_z', v_\theta\) and \(v_z'\) are smaller than <tt>1e-3</tt> on-axis.
  !! Nothing is checked if the geometry is Cartesian.
  !! @warning Throws a warning if:
  !!
  !! - \(B_\theta\) is not zero on-axis.
  !! - \(B_z'\) is not zero on-axis.
  !! - \(v_\theta\) is not zero on-axis.
  !! - \(v_z'\) is not zero on-axis. @endwarning
  subroutine check_on_axis_values(settings, background)
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background

    if (settings%grid%get_geometry() == "Cartesian") return
    if (settings%grid%get_grid_start() > 0.0_dp) return

    ! LCOV_EXCL_START
    if (.not. is_zero(background%magnetic%B02(0.0_dp))) then
      call logger%warning( &
        "B_theta non-zero on axis! Value: " // str(background%magnetic%B02(0.0_dp)) &
      )
    end if
    if (.not. is_zero(background%magnetic%dB03(0.0_dp))) then
      call logger%warning( &
        "dBz/dr non-zero on axis! Value: " // str(background%magnetic%dB03(0.0_dp)) &
      )
    end if
    if (.not. is_zero(background%velocity%v02(0.0_dp))) then
      call logger%warning( &
        "v_theta non-zero on axis! Value: " // str(background%velocity%v02(0.0_dp)) &
      )
    end if
    if (.not. is_zero(background%velocity%dv03(0.0_dp))) then
      call logger%warning( &
        "dvz_dr non-zero on axis! Value: " // str(background%velocity%dv03(0.0_dp)) &
      )
    end if
    ! LCOV_EXCL_STOP
  end subroutine check_on_axis_values


  !> Checks the standard force-balance equation for the equilibrium state. This
  !! results in three expressions,
  !! $$ \Bigl(p_0 + \frac{1}{2}B_0^2\Bigr)' + \rho_0 g + \rho_0 v_{01} v_{01}'
  !!    - \frac{\varepsilon'}{\varepsilon}\bigl(\rho_0 v_{02}^2 - B_{02}^2\bigr) = 0, $$
  !! $$ \rho_0 v_{01} \bigl( v_{02}' + \frac{\varepsilon'}{\varepsilon} v_{02} \bigr)
  !!    - \frac{B_{01}}{\varepsilon} (\varepsilon B_{02})' = 0, $$
  !! $$ \rho_0 v_{01} v_{03}' - B_{01} B_{03}' = 0, $$
  !! and they should all be fulfilled.
  !! @warning   Throws a warning if force-balance is not satisfied.
  subroutine standard_equil_conditions(settings, background, grav_field)
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    !> the type containing the gravity attributes
    type(gravity_type), intent(in)      :: grav_field

    real(dp) :: x
    real(dp) :: rho0, drho0, B01, B02, dB02, B03, dB03, T0, dT0, grav
    real(dp) :: v01, v02, v03, dv01, dv02, dv03
    real(dp) :: eps, d_eps, r(3), discrepancy(3)
    real(dp) :: eq_cond(settings%grid%get_gauss_gridpts(), 3)
    integer :: i, j, counter(3)
    logical :: satisfied(3)

    satisfied = .true.
    discrepancy = 0.0d0
    counter = 0
    do i = 1, settings%grid%get_gauss_gridpts()
      x = grid_gauss(i)
      rho0 = background%density%rho0(x)
      drho0 = background%density%drho0(x)
      B01 = background%magnetic%B01(x)
      B02 = background%magnetic%B02(x)
      B03 = background%magnetic%B03(x)
      dB02 = background%magnetic%dB02(x)
      dB03 = background%magnetic%dB03(x)
      T0 = background%temperature%T0(x)
      dT0 = background%temperature%dT0(x)
      v01 = background%velocity%v01(x)
      v02 = background%velocity%v02(x)
      v03 = background%velocity%v03(x)
      dv01 = background%velocity%dv01(x)
      dv02 = background%velocity%dv02(x)
      dv03 = background%velocity%dv03(x)
      grav = grav_field % grav(i)
      eps = eps_grid(i)
      d_eps = d_eps_grid_dr(i)

      eq_cond(i, 1) = ( &
        drho0 * T0 &
        + rho0 * dT0 &
        + B02 * dB02 &
        + B03 * dB03 &
        + rho0 * grav &
        - (d_eps/eps) * (rho0 * v02**2 - B02**2) &
        + rho0 * v01 * dv01 &
      )
      eq_cond(i, 2) = ( &
        rho0 * v01 * (dv02 + v02 * d_eps / eps) - B01 * (dB02 + B02 * d_eps / eps) &
      )
      eq_cond(i, 3) = rho0 * v01 * dv03 - B01 * dB03

      do j = 1, 3
        if (abs(eq_cond(i, j)) > dp_LIMIT) then
          counter(j) = counter(j) + 1
          satisfied(j) = .false.
          if (abs(eq_cond(i, j)) > discrepancy(j)) then
            discrepancy(j) = abs(eq_cond(i, j))
            r(j) = grid_gauss(i)
          end if
        end if
      end do
    end do

    do j = 1, 3
      if (satisfied(j)) then
        cycle
      end if
      ! LCOV_EXCL_START
      call logger%warning("standard equilibrium conditions not satisfied!")
      call logger%warning( &
        "location of largest discrepancy (" // str(j) // "): x = " // str(r(j)) &
      )
      call logger%warning( &
        "value of largest discrepancy (" // str(j) // "): " &
        // str(discrepancy(j), fmt=exp_fmt) &
      )
      call logger%warning( &
        "amount of nodes not satisfying criterion (" // str(j) // "): " &
        // str(counter(j)) &
      )
      ! LCOV_EXCL_STOP
    end do
  end subroutine standard_equil_conditions


  !> Enforces the non-adiabatic force-balance equation for the equilibrium state.
  !! This is given by
  !! $$
  !! T_0 \rho_0 \frac{\left(\varepsilon v_{01}\right)'}{\varepsilon}
  !! + \rho_0 \mathscr{L}_0
  !! - B_{01}^2\left[\frac{\kappa_{\parallel,0} - \kappa_{\perp,0}}{B_0^2} T_0'\right]'
  !! - \frac{1}{\varepsilon}\left(\varepsilon \kappa_{\perp, 0} T_0'\right)'
  !! + \frac{1}{(\gamma - 1)}T_0'\rho_0 v_{01} = 0
  !! $$
  subroutine check_energy_balance(settings, background, rc_field, kappa_field)
    !> the settings object
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    !> the type containing the radiative cooling attributes
    type(cooling_type), intent(inout)  :: rc_field
    !> the type containing the thermal conduction attributes
    type(conduction_type), intent(in) :: kappa_field

    real(dp) :: x, rho0, drho0, T0, dT0, ddT0
    real(dp) :: B01, B02, dB02, B03, dB03, B0
    real(dp) :: v01, dv01
    real(dp) :: kappa_perp, dkappa_perp_dr, Kp, dKp
    real(dp) :: eps, deps
    real(dp) :: balance_terms
    real(dp) :: discrepancy, r, eq_cond
    logical :: satisfied
    integer :: i, counter

    satisfied = .true.
    discrepancy = 0.0_dp
    counter = 0
    do i = 1, settings%grid%get_gauss_gridpts()
      x = grid_gauss(i)
      rho0 = background%density%rho0(x)
      drho0 = background%density%drho0(x)
      B01 = background%magnetic%B01(x)
      B02 = background%magnetic%B02(x)
      B03 = background%magnetic%B03(x)
      dB02 = background%magnetic%dB02(x)
      dB03 = background%magnetic%dB03(x)
      B0 = background%magnetic%get_B0(x)
      T0 = background%temperature%T0(x)
      dT0 = background%temperature%dT0(x)
      ddT0 = background%temperature%ddT0(x)
      v01 = background%velocity%v01(x)
      dv01 = background%velocity%dv01(x)
      eps = eps_grid(i)
      deps = d_eps_grid_dr(i)
      kappa_perp = kappa_field % kappa_perp(i)
      dkappa_perp_dr = kappa_field % d_kappa_perp_dr(i)
      Kp = kappa_field % prefactor(i)
      dKp = kappa_field % d_prefactor_dr(i)

      balance_terms = ( &
        -(T0 / eps) * (deps * v01 + eps * dv01) &
        + dT0 * v01 / settings%physics%get_gamma_1() &
        + (B01**2 / rho0) * (dKp * dT0 + Kp * ddT0) &
        + ( &
          deps * kappa_perp * dT0 &
          + eps * dkappa_perp_dr * dT0 &
          + eps * kappa_perp * ddT0 &
        ) / (eps * rho0) &
      )

      if (settings%physics%cooling%ensure_equilibrium) then
        rc_field%L0(i) = balance_terms
      else
        rc_field%L0(i) = rho0 * rc_field%lambdaT(i) - rc_field%H0(i)
      end if

      eq_cond = ( &
        T0 * rho0 * (deps * v01 + eps * dv01) / eps &
        + rho0 * rc_field%L0(i) &
        - B01**2 * (Kp * dT0 + dKp * T0) &
        - (1.0_dp / eps) * ( &
          deps * kappa_perp * dT0 &
          + eps * dkappa_perp_dr * dT0 &
          + eps * kappa_perp * ddT0 &
        ) &
        + (1.0_dp / settings%physics%get_gamma_1()) * dT0 * rho0 * v01 &
      )
      if (.not. is_zero(eq_cond)) then
        counter = counter + 1
        satisfied = .false.
        if (abs(eq_cond) > discrepancy) then
          discrepancy = abs(eq_cond)
          r = grid_gauss(i)
        end if
      end if
    end do

    ! LCOV_EXCL_START
    if (.not. satisfied) then
      call logger%warning("energy balance not satisfied!")
      call logger%warning("location of largest discrepancy: x = " // str(r))
      call logger%warning( &
        "value of largest discrepancy: " // str(discrepancy, fmt=exp_fmt) &
      )
      call logger%warning("amount of nodes not satisfying criterion: " // str(counter))
    end if
    ! LCOV_EXCL_STOP
  end subroutine check_energy_balance


  !> Checks the induction equation for the equilibrium state. The two (nonzero)
  !! resulting expressions are
  !! $$ (B_{01} v_{02} - (B_{02} v_{01})' = 0, $$
  !! $$
  !! \frac{1}{\varepsilon} \bigl(\varepsilon (B_{01}v_{03} - B_{03}v_{01}) \bigr)' = 0
  !! $$
  !! and should both be fulfilled.
  !! @warning   Throws a warning if the equilibrium induction equation is not satisfied.
  subroutine induction_equil_conditions(settings, background)
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background

    real(dp)  :: x, B01, B02, dB02, B03, dB03, v01, v02, v03, dv01, dv02, dv03
    real(dp)  :: eps, d_eps, r(2), discrepancy(2)
    real(dp)  :: eq_cond(settings%grid%get_gauss_gridpts(), 2)
    integer   :: i, j, counter(2)
    logical   :: satisfied(2)

    satisfied = .true.
    discrepancy = 0.0d0
    counter = 0
    do i = 1, settings%grid%get_gauss_gridpts()
      x = grid_gauss(i)
      B01 = background%magnetic%B01(x)
      B02 = background%magnetic%B02(x)
      B03 = background%magnetic%B03(x)
      dB02 = background%magnetic%dB02(x)
      dB03 = background%magnetic%dB03(x)
      v01 = background%velocity%v01(x)
      v02 = background%velocity%v02(x)
      v03 = background%velocity%v03(x)
      dv01 = background%velocity%dv01(x)
      dv02 = background%velocity%dv02(x)
      dv03 = background%velocity%dv03(x)
      eps = eps_grid(i)
      d_eps = d_eps_grid_dr(i)

      eq_cond(i, 1) = B01 * dv02 - B02 * dv01 - dB02 * v01
      eq_cond(i, 2) = ( &
        B01 * dv03 - dB03 * v01 - B03 * dv01 + d_eps * (B01 * v03 - B03 * v01) / eps &
      )
      do j = 1, 2
        if (abs(eq_cond(i, j)) > dp_LIMIT) then
          counter(j) = counter(j) + 1
          satisfied(j) = .false.
          if (abs(eq_cond(i, j)) > discrepancy(j)) then
            discrepancy(j) = abs(eq_cond(i, j))
            r(j) = grid_gauss(i)
          end if
        end if
      end do
    end do

    do j = 1, 2
      if (satisfied(j)) then
        cycle
      end if
      ! LCOV_EXCL_START
      call logger%warning("induction equilibrium conditions not satisfied!")
      call logger%warning( &
        "location of largest discrepancy (" // str(j) // "): x = " // str(r(j)) &
      )
      call logger%warning( &
        "value of largest discrepancy (" // str(j) // "): " &
        // str(discrepancy(j), fmt=exp_fmt) &
      )
      call logger%warning( &
        "amount of nodes not satisfying criterion (" // str(j) // "): " &
        // str(counter(j)) &
      )
    end do
    ! LCOV_EXCL_STOP
  end subroutine induction_equil_conditions



  !> Checks the continuity equation for the equilibrium state. This is given by
  !! $$ \frac{1}{\varepsilon} \bigl( \varepsilon \rho_0 v_{01} \bigr)' = 0. $$
  !! @warning   Throws a warning if equilibrium continuity is not satisfied.
  subroutine continuity_equil_conditions(settings, background)
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background

    real(dp)  :: x, rho0, drho0, v01, dv01
    real(dp)  :: eps, d_eps, r, discrepancy
    real(dp)  :: eq_cond(settings%grid%get_gauss_gridpts())
    integer   :: i, counter
    logical   :: satisfied

    satisfied = .true.
    discrepancy = 0.0d0
    counter = 0
    do i = 1, settings%grid%get_gauss_gridpts() - 1
      x = grid_gauss(i)
      rho0 = background%density%rho0(x)
      drho0 = background%density%drho0(x)
      v01 = background%velocity%v01(x)
      dv01 = background%velocity%dv01(x)
      eps = eps_grid(i)
      d_eps = d_eps_grid_dr(i)

      eq_cond(i) = drho0 * v01 + rho0 * dv01 + rho0 * v01 * d_eps / eps

      if (abs(eq_cond(i)) > dp_LIMIT) then
        counter = counter + 1
        satisfied = .false.
        if (abs(eq_cond(i)) > discrepancy) then
          discrepancy = abs(eq_cond(i))
          r = grid_gauss(i)
        end if
      end if
    end do

    ! LCOV_EXCL_START
    if (.not. satisfied) then
      call logger%warning("continuity equilibrium conditions not satisfied!")
      call logger%warning("location of largest discrepancy: x = " // str(r))
      call logger%warning( &
        "value of largest discrepancy: " // str(discrepancy, fmt=exp_fmt) &
      )
      call logger%warning("amount of nodes not satisfying criterion: " // str(counter))
    end if
    ! LCOV_EXCL_STOP
  end subroutine continuity_equil_conditions

end module mod_inspections
