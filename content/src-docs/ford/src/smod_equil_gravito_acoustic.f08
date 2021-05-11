! =============================================================================
!> This submodule defines an equilibrium in Cartesian geometry with a
!! stratified equilibrium profile, giving rise to gravito-acoustic waves.
!! No magnetic fields are included, such that this treats the hydrodynamic regime.
!! The geometry can be overridden using the parfile.
!!
!! This equilibrium is taken from section 7.2.3, p. 242 in
!! _Goedbloed, H., Keppens, R., & Poedts, S. (2019). Magnetohydrodynamics of Laboratory
!!  and Astrophysical Plasmas. Cambridge University Press._ [DOI](http://doi.org/10.1017/9781316403679).
!! @note Default values are given by
!!
!! - <tt>k2</tt> = \(\pi\)
!! - <tt>k3</tt> = \(\pi\)
!! - <tt>cte_p0</tt> = 1 : used to set the pressure value.
!! - <tt>alpha</tt> = 20.42 : used to constrain the density.
!! - <tt>g</tt> = 0.5 : used to set the gravity constant.
!!
!! and can all be changed in the parfile. @endnote
!! @warning This equilibrium has no regression test yet! @endwarning
submodule (mod_equilibrium) smod_equil_gravito_acoustic
  implicit none

contains

  module subroutine gravito_acoustic_eq()
    use mod_equilibrium_params, only: g, cte_rho0, cte_p0, alpha

    real(dp)  :: x, g
    integer   :: i

    call allow_geometry_override(default_geometry='Cartesian', default_x_start=0.0d0, default_x_end=1.0d0)
    call initialise_grid()

    if (use_defaults) then
      external_gravity = .true.

      k2 = dpi
      k3 = dpi
      cte_p0 = 1.0d0
      alpha = 20.42d0
      g = 0.5d0
    end if

    cte_rho0 = alpha * cte_p0 / g
    T_field % T0      = cte_p0 / cte_rho0
    grav_field % grav = g

    do i = 1, gauss_gridpts
      x = grid_gauss(i)

      rho_field % rho0(i) = cte_rho0 * exp(-alpha*x)
      rho_field % d_rho0_dr(i) = -alpha * (rho_field % rho0(i))
    end do

  end subroutine gravito_acoustic_eq

end submodule smod_equil_gravito_acoustic
