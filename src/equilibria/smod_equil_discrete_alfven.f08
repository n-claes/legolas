! =============================================================================
!> This submodule defines an equilibrium in cylindrical geometry with
!! an axial current profile (\(\nu = 2\)), modelling a solar coronal loop
!! in which discrete Alfvén waves are present. The geometry can be overridden in the parfile.
!!
!! This equilibrium is taken from
!! _Keppens, Rony, Ronald AM Van Der Linden, and Marcel Goossens.
!! "Non-adiabatic discrete Alfven waves in coronal loops and prominences.",
!! Solar physics 144.2 (1993): 267-281_.
!!
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 1
!! - <tt>k3</tt> = 0.05
!! - <tt>j0</tt> = 0.125 : used to set the current.
!! - <tt>delta</tt> = 0.2 : used in the density profile.
!! - cooling_curve = 'rosner'
!! - parallel thermal conduction, no perpendicular conduction
!!
!! and normalisations given by
!!
!! - <tt>unit_density</tt> = 1.5e-15 gcm-3
!! - <tt>unit_magneticfield</tt> = 50 Gauss
!! - <tt>unit_length</tt> = 1e10 cm
!!
!! and can all be changed in the parfile. @endnote
submodule (mod_equilibrium) smod_equil_discrete_alfven
  implicit none

contains

  module procedure discrete_alfven_eq
    use mod_equilibrium_params, only: j0, delta

    real(dp) :: r, x_end
    real(dp), allocatable :: p_r(:), dp_r(:)
    integer :: i, gauss_gridpts

    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%grid%set_geometry("cylindrical")
      call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)
      call settings%physics%enable_cooling(cooling_curve="rosner")
      call settings%physics%enable_parallel_conduction()
      call settings%units%set_units_from_density( &
        unit_density=1.5d-15, &
        unit_magneticfield=50.0d0, &
        unit_length=1.0d10, &
        mean_molecular_weight=1.0d0 & ! pure proton plasma
      )

      j0 = 0.125d0
      delta = 0.2d0   ! d parameter in density prescription
      k2 = 1.0d0
      k3 = 0.05d0
    end if ! LCOV_EXCL_STOP
    gauss_gridpts = settings%grid%get_gauss_gridpts()
    allocate(p_r(gauss_gridpts))
    allocate(dp_r(gauss_gridpts))

    call initialise_grid(settings)
    x_end = settings%grid%get_grid_end()

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      rho_field % rho0(i) = 1.0d0 - (1.0d0 - delta) * (r / x_end)**2
      B_field % B02(i) = j0 * r * (r**4 - 3.0d0 * r**2 + 3.0d0) / 6.0d0
      B_field % B03(i) = 1.0d0
      B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      p_r(i) = (j0**2 / 720.0d0) * ( &
        12.0d0 * (x_end**10 - r**10) &
        - 75.0d0 * (x_end**8 - r**8) &
        + 200.0d0 * (x_end**6 - r**6) &
        - 270.0d0 * (x_end**4 - r**4) &
        + 180.0d0 * (x_end**2 - r**2) &
      )
      T_field % T0(i) = p_r(i) / (rho_field % rho0(i))

      rho_field % d_rho0_dr(i) = 2.0d0 * r * (delta - 1.0d0) / x_end**2
      B_field % d_B02_dr(i) = j0 * (5.0d0 * r**4 - 9.0d0 * r**2 + 3.0d0) / 6.0d0
      dp_r(i) = j0**2 * r * ( &
        -r**8 + 5.0d0 * r**6 - 10.0d0 * r**4 + 9.0d0 * r**2 - 3.0d0 &
      ) / 6.0d0
      T_field % d_T0_dr(i) = ( &
        dp_r(i) * (rho_field % rho0(i)) - (rho_field % d_rho0_dr(i)) * p_r(i) &
      ) / (rho_field % rho0(i))**2
    end do
    deallocate(p_r)
    deallocate(dp_r)
  end procedure discrete_alfven_eq

end submodule smod_equil_discrete_alfven
