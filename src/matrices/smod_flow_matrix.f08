submodule (mod_matrix_manager) smod_flow_matrix
  implicit none

contains

  module procedure add_flow_matrix_terms
    real(dp)  :: eps, deps
    real(dp)  :: rho, drho
    real(dp)  :: T0
    real(dp)  :: v01, dv01, drv01
    real(dp)  :: v02, dv02, drv02
    real(dp)  :: v03, dv03
    real(dp)  :: Vop
    real(dp) :: gamma_1
    type(matrix_elements_t) :: elements

    gamma_1 = settings%physics%get_gamma_1()
    ! grid variables
    eps = grid%get_eps(x_gauss)
    deps = grid%get_deps()
    ! density variables
    rho = background%density%rho0(x_gauss)
    drho = background%density%drho0(x_gauss)
    ! temperature variables
    T0 = background%temperature%T0(x_gauss)
    ! flow variables
    v01 = background%velocity%v01(x_gauss)
    dv01 = background%velocity%dv01(x_gauss)
    drv01 = deps * v01 + eps * dv01
    v02 = background%velocity%v02(x_gauss)
    dv02 = background%velocity%dv02(x_gauss)
    drv02 = deps * v02 + eps * dv02
    v03 = background%velocity%v03(x_gauss)
    dv03 = background%velocity%dv03(x_gauss)
    Vop = k2 * v02 / eps + k3 * v03

    elements = new_matrix_elements(state_vector=settings%get_state_vector())

    ! ==================== Quadratic * Quadratic ====================
    call elements%add(Vop - ic * dv01, "rho", "rho", spline1=h_quad, spline2=h_quad)
    call elements%add( &
      -drv02 * ic * v01 / eps, "v2", "rho", spline1=h_quad, spline2=h_quad &
    )
    call elements%add( &
      rho * (eps * Vop - ic * deps * v01), "v2", "v2", spline1=h_quad, spline2=h_quad &
    )
    call elements%add(-ic * v01 * dv03, "v3", "rho", spline1=h_quad, spline2=h_quad)
    call elements%add( &
      rho * (Vop + ic * dv01) + (deps * rho / eps + drho) * ic * v01, &
      "v3", &
      "v3", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    call elements%add(eps * Vop, "a1", "a1", spline1=h_quad, spline2=h_quad)

    ! ==================== Quadratic * dQuadratic ====================
    call elements%add(-ic * v01, "rho", "rho", spline1=h_quad, spline2=dh_quad)
    call elements%add( &
      -ic * eps * rho * v01, "v2", "v2", spline1=h_quad, spline2=dh_quad &
    )

    ! ==================== Cubic * Quadratic ====================
    call elements%add( &
      v01 * dv01 - deps * v02**2 / eps, "v1", "rho", spline1=h_cubic, spline2=h_quad &
    )
    call elements%add( &
      -2.0d0 * deps * rho * v02, "v1", "v2", spline1=h_cubic, spline2=h_quad &
    )
    call elements%add(ic * v01 * k2, "a2", "a1", spline1=h_cubic, spline2=h_quad)
    call elements%add(eps * k3 * ic * v01, "a3", "a1", spline1=h_cubic, spline2=h_quad)

    ! ==================== Cubic * Cubic ====================
    call elements%add( &
      rho * Vop + (deps * rho / eps + drho) * ic * v01, &
      "v1", &
      "v1", &
      spline1=h_cubic, &
      spline2=h_cubic &
    )
    call elements%add(k3 * v03, "a2", "a2", spline1=h_cubic, spline2=h_cubic)
    call elements%add(-k2 * v03, "a2", "a3", spline1=h_cubic, spline2=h_cubic)
    call elements%add(-k3 * v02, "a3", "a2", spline1=h_cubic, spline2=h_cubic)
    call elements%add(k2 * v02, "a3", "a3", spline1=h_cubic, spline2=h_cubic)

    ! ==================== dCubic * Cubic ====================
    call elements%add(ic * rho * v01, "v1", "v1", spline1=dh_cubic, spline2=h_cubic)

    ! ==================== Quadratic * Cubic ====================
    call elements%add(-drv02 * rho / eps, "v2", "v1", spline1=h_quad, spline2=h_cubic)
    call elements%add(-rho * dv03, "v3", "v1", spline1=h_quad, spline2=h_cubic)

    ! ==================== dQuadratic * Quadratic ====================
    call elements%add(ic * rho * v01, "v3", "v3", spline1=dh_quad, spline2=h_quad)

    ! ==================== Quadratic * dCubic ====================
    call elements%add(-v02, "a1", "a2", spline1=h_quad, spline2=dh_cubic)
    call elements%add(-eps * v03, "a1", "a3", spline1=h_quad, spline2=dh_cubic)

    ! ==================== Cubic * dCubic ====================
    call elements%add(-ic * v01, "a2", "a2", spline1=h_cubic, spline2=dh_cubic)
    call elements%add(-ic * v01, "a3", "a3", spline1=h_cubic, spline2=dh_cubic)

    if (.not. settings%physics%is_incompressible) then
      ! ==================== Quadratic * Quadratic ====================
      call elements%add( &
        -ic * gamma_1 * drv01 * T0 / eps, "T", "rho", spline1=h_quad, spline2=h_quad &
      )
      call elements%add( &
        rho * (Vop + ic * dv01 - ic * gamma_1 * drv01 / eps) &
          + ic * v01 * (deps * rho / eps + drho), &
        "T", &
        "T", &
        spline1=h_quad, &
        spline2=h_quad &
      )
      ! ==================== dQuadratic * Quadratic ====================
      call elements%add(ic * rho * v01, "T", "T", spline1=dh_quad, spline2=h_quad)
    end if

    call add_to_quadblock(quadblock, elements, weight, settings%dims)
    call elements%delete()
  end procedure add_flow_matrix_terms

end submodule smod_flow_matrix
