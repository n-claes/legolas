! =============================================================================
!> This submodule defines a magnetic flux tube embedded in a uniform magnetic
!! environment. In this case the flux tube is under coronal conditions
!! \( c_s < c_A < c_{Ae} \) where the subscript e denotes the outer region.
!! More specifically the equilibrium is defined as
!! \( c_{Ae} = 5c_s, c_{se} = c_s/2, c_A = 2c_s \). The geometry can be overridden
!! in the parfile, and is cylindrical by default for \( r \in [0, 10] \).
!!
!! This equilibrium is taken from chapter 6, fig. 6.7 in
!! _Roberts, Bernard (2019). MHD Waves in the Solar Atmosphere.
!! Cambridge University Press._ [DOI](https://doi.org/10.1017/9781108613774).
!! @note For best results, it is recommended to enable mesh accumulation. @endnote
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 0
!! - <tt>k3</tt> = 2
!! - <tt>cte_rho0</tt> = 1 : density value for the inner tube.
!! - <tt>cte_p0</tt> = 1 : pressure value for the inner tube.
!! - <tt>r0</tt> = 1 : radius of the inner tube.
!!
!! and can all be changed in the parfile. @endnote
! SUBMODULE: smod_equil_coronal_flux_tube
submodule(mod_equilibrium) smod_equil_coronal_flux_tube
  implicit none

contains

  module procedure coronal_flux_tube_eq
    use mod_global_variables, only: dp_LIMIT
    use mod_equilibrium_params, only: cte_rho0, cte_p0, r0

    real(dp)  :: r, rho_e, p_e, B_0, B_e
    real(dp), allocatable :: custom_grid(:)
    real(dp)  :: width, a_l, a_r, pct1, pct2, pct3, dx, dx1, dx2, dx3
    real(dp) :: gamma
    real(dp) :: x_start, x_end
    integer :: gridpts
    integer :: i, N1, N2, N3

    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%grid%set_geometry("cylindrical")
      call settings%grid%set_grid_boundaries(0.0_dp, 10.0_dp)

      cte_rho0 = 1.0d0
      cte_p0 = 1.0d0
      r0 = 1.0d0

      k2 = 0.0d0
      k3 = 2.0d0
    end if ! LCOV_EXCL_STOP
    gamma = settings%physics%get_gamma()
    gridpts = settings%grid%get_gridpts()
    x_start = settings%grid%get_grid_start()
    x_end = settings%grid%get_grid_end()

    ! width of transition region
    width = 0.1d0
    ! start of transition region
    a_l = r0 - width / 2.0d0
    ! end of transition region
    a_r = a_l + width
    pct1 = 0.6d0  ! percentage of points in inner tube
    pct2 = 0.3d0  ! percentage of points in transition region
    pct3 = 0.1d0  ! percentage of points in outer tube
    ! for a custom grid we have to manually avoid setting r0 = 0
    if ( &
      settings%grid%get_geometry() == "cylindrical" &
      .and. .not. settings%grid%force_r0 &
    ) then
      x_start = 2.5d-2
    end if
    N1 = int(pct1 * gridpts)
    dx1 = (a_l - x_start) / (N1 - 1)  ! -1 since first point is x0
    N2 = int(pct2 * gridpts)
    dx2 = (a_r - a_l) / N2
    N3 = int(pct3 * gridpts)
    dx3 = (x_end - a_r) / N3
    ! fill the grid
    allocate(custom_grid(gridpts))
    custom_grid(1) = x_start
    do i = 2, gridpts
      if (i <= N1) then
        dx = dx1
      else if (N1 < i .and. i <= N1 + N2) then
        dx = dx2
      else
        dx = dx3
      end if
      custom_grid(i) = custom_grid(i - 1) + dx
    end do
    call initialise_grid(settings, custom_grid=custom_grid)
    deallocate(custom_grid)

    rho_e = 4.0d0 * (2.0d0 * gamma + 1.0d0) * cte_rho0 / (50.0d0 * gamma + 1.0d0)
    p_e = (2.0d0 * gamma + 1.0d0) * cte_p0 / (50.0d0 * gamma + 1.0d0)
    B_0 = 2.0d0 * sqrt(gamma * cte_p0)
    B_e = 10.0d0 * sqrt( &
      gamma * cte_p0 * (2.0d0 * gamma + 1.0d0) / (50.0d0 * gamma + 1.0d0) &
    )

    if (r0 > x_end) then
      call log_message("equilibrium: inner cylinder radius r0 > x_end", level="error")
    else if (r0 < x_start) then
      call log_message("equilibrium: inner cylinder radius r0 < x_start", level="error")
    end if

    ! check pressure balance
    if (abs(cte_p0 + 0.5d0 * B_0**2 - p_e - 0.5d0 * B_e**2) > dp_LIMIT) then
      call log_message( &
        "equilibrium: total pressure balance not satisfied", level="error" &
      )
    end if

    do i = 1, settings%grid%get_gauss_gridpts()
      r = grid_gauss(i)

      if (r > r0) then
        rho_field % rho0(i) = rho_e
        B_field % B02(i) = 0.0d0
        B_field % B03(i) = B_e
        T_field % T0(i) = p_e / rho_e
      else
        rho_field % rho0(i) = cte_rho0
        B_field % B02(i) = 0.0d0
        B_field % B03(i) = B_0
        T_field % T0(i) = cte_p0 / cte_rho0
      end if
      B_field % B0(i)  = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
    end do

  end procedure coronal_flux_tube_eq

end submodule smod_equil_coronal_flux_tube
