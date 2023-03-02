submodule (mod_matrix_manager) smod_flow_matrix
  use mod_equilibrium, only: v_field
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
    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)
    ! density variables
    rho = rho_field % rho0(gauss_idx)
    drho = rho_field % d_rho0_dr(gauss_idx)
    ! temperature variables
    T0 = T_field % T0(gauss_idx)
    ! flow variables
    v01 = v_field % v01(gauss_idx)
    dv01 = v_field % d_v01_dr(gauss_idx)
    drv01 = deps * v01 + eps * dv01
    v02 = v_field % v02(gauss_idx)
    dv02 = v_field % d_v02_dr(gauss_idx)
    drv02 = deps * v02 + eps * dv02
    v03 = v_field % v03(gauss_idx)
    dv03 = v_field % d_v03_dr(gauss_idx)
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

    call add_to_quadblock(quadblock, elements, current_weight, settings%dims)
    call elements%delete()
  end procedure add_flow_matrix_terms

end submodule smod_flow_matrix
