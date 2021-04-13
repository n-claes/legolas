! =============================================================================
!> Assembly of the finite element matrices \(\mathcal{A}\) and \(\mathcal{B}\) in the matrix
!! eigenvalue problem \(\mathcal{A}\textbf{X} = \omega\mathcal{B}\textbf{X}\).
!! Matrix creation proceeds by first assembling a quadblock at every grid interval
!! using the integrals and finite element basis functions, after which this quadblock
!! is shifted along the main diagonal.
!! @note    Boundary conditions are imposed after matrix assembly is completed.
module mod_matrix_creation
  use mod_global_variables, only: dp, gridpts, matrix_gridpts, dim_quadblock, &
                                      dim_subblock, thermal_conduction
  implicit none

  private

  !> array containing the 4 quadratic basis functions
  real(dp)  :: h_quadratic(4)
  !> array containing the 4 cubic basis functions
  real(dp)  :: h_cubic(4)
  !> array containing the derivatives of the 4 quadratic basis functions
  real(dp)  :: dh_quadratic_dr(4)
  !> array containing the derivatives of the 4 cubic basis functions
  real(dp)  :: dh_cubic_dr(4)

  public :: create_matrices

contains


  !> Main subroutine to assemble the matrices A and B, which are
  !! already allocated when entering this subroutine but not yet initialised.
  !! The quadblock is calculated for every grid interval and used to assemble
  !! both matrices. On exit, both matrices are fully assembled and
  !! boundary conditions are imposed.
  subroutine create_matrices(matrix_B, matrix_A)
    use mod_global_variables, only: gaussian_weights, n_gauss, &
                                    ic, viscosity, viscosity_value, &
                                    hall_mhd, elec_inertia
    use mod_grid, only: grid, grid_gauss, eps_grid, d_eps_grid_dr
    use mod_spline_functions, only: quadratic_factors, quadratic_factors_deriv, &
                                    cubic_factors, cubic_factors_deriv
    use mod_boundary_conditions, only: apply_boundary_conditions
    use mod_viscosity, only: get_viscosity_terms
    use mod_hallmhd, only: get_hallterm, get_hallterm_Bmat
    use mod_equilibrium, only: hall_field

    !> the B-matrix
    real(dp), intent(inout)     :: matrix_B(matrix_gridpts, matrix_gridpts)
    !> the A-matrix
    complex(dp), intent(inout)  :: matrix_A(matrix_gridpts, matrix_gridpts)
    complex(dp)                 :: quadblock_B(dim_quadblock, dim_quadblock)
    complex(dp)                 :: quadblock_A(dim_quadblock, dim_quadblock)
    complex(dp)                 :: quadblock_viscosity(dim_quadblock, dim_quadblock)
    complex(dp)                 :: quadblock_Hall(dim_quadblock, dim_quadblock)
    complex(dp)                 :: quadblock_HallB(dim_quadblock, dim_quadblock)

    real(dp)                    :: r, r_lo, r_hi, eps, d_eps_dr, curr_weight
    integer                     :: i, j, gauss_idx, k, l
    integer                     :: quadblock_idx, idx1, idx2
    real(dp)                    :: hf, inf

    ! initialise matrices (A is complex, B is real)
    matrix_B = 0.0d0
    matrix_A = (0.0d0, 0.0d0)
    ! initialise quadblock index to 0 (used to shift the block along diagonal)
    quadblock_idx = 0

    do i = 1, gridpts-1
      ! reset quadblock (complex)
      quadblock_B = (0.0d0, 0.0d0)
      quadblock_A = (0.0d0, 0.0d0)

      ! This integrates in the interval grid(i) to grid(i + 1)
      r_lo = grid(i)
      r_hi = grid(i+1)

      ! in every grid interval loop over gaussian points to calculate integral
      do j = 1, n_gauss
        ! current grid index (in Gaussian grid)
        gauss_idx = (i - 1)*n_gauss + j
        r = grid_gauss(gauss_idx)

        ! define scale factor
        eps = eps_grid(gauss_idx)
        d_eps_dr = d_eps_grid_dr(gauss_idx)

        curr_weight = gaussian_weights(j)

        ! calculate spline functions for this point (r) in the grid interval
        call quadratic_factors(r, r_lo, r_hi, h_quadratic)
        call quadratic_factors_deriv(r, r_lo, r_hi, dh_quadratic_dr)
        call cubic_factors(r, r_lo, r_hi, h_cubic)
        call cubic_factors_deriv(r, r_lo, r_hi, dh_cubic_dr)

        ! calculate matrix elements for this point
        call get_B_elements(gauss_idx, eps, curr_weight, quadblock_B)
        call get_A_elements(gauss_idx, eps, d_eps_dr, curr_weight, quadblock_A)

        if (viscosity) then
          call get_viscosity_terms(gauss_idx, eps, d_eps_dr, curr_weight, quadblock_viscosity, &
                                    h_quadratic, dh_quadratic_dr, h_cubic, dh_cubic_dr)
          quadblock_A = quadblock_A + ic * viscosity_value * quadblock_viscosity
        end if
        if (hall_mhd) then
          call get_hallterm(gauss_idx, eps, d_eps_dr, curr_weight, quadblock_Hall, &
                            h_quadratic, dh_quadratic_dr, h_cubic, dh_cubic_dr)
          hf = hall_field % hallfactor(gauss_idx)
          quadblock_A = quadblock_A - hf * quadblock_Hall
        end if

        if (elec_inertia) then
          call get_hallterm_Bmat(gauss_idx, eps, d_eps_dr, curr_weight, quadblock_HallB, &
                            h_quadratic, h_cubic, dh_cubic_dr)
          inf = hall_field % inertiafactor(gauss_idx)
          quadblock_B = quadblock_B + inf * quadblock_HallB
        end if
      end do   ! ends iteration gaussian points

      ! multiply by dx for integral
      quadblock_B = quadblock_B * (r_hi - r_lo)
      quadblock_A = quadblock_A * (r_hi - r_lo)

      ! fill matrices with quadblock entries.
      do k = 1, dim_quadblock
        do l = 1, dim_quadblock
          ! quadblock is shifted on main diagonal using idx1 and idx2
          ! dimension dim_subblock is added instead of dim_quadblock, as the
          ! bottom-right part of the quadblock overlaps with the top-left part
          ! when shifting along the main diagonal
          idx1 = k + quadblock_idx
          idx2 = l + quadblock_idx
          matrix_B(idx1, idx2) = matrix_B(idx1, idx2) + real(quadblock_B(k, l))
          matrix_A(idx1, idx2) = matrix_A(idx1, idx2) + quadblock_A(k, l)
        end do
      end do
      quadblock_idx = quadblock_idx + dim_subblock
    end do    ! ends iteration over grid points

    ! ! apply boundary conditions on assembled matrices
    call apply_boundary_conditions(matrix_A, matrix_B)
  end subroutine create_matrices


  !> Retrieves the B-matrix elements for a given quadblock corresponding
  !! to a particular point in the grid and a particular Gaussian weight.
  !! This routine is called <tt>n_gauss</tt> times for every grid interval.
  !! @note    When this subroutine gets called the basis functions defined at module scope
  !!          have already been recalculated and correspond to this particular value
  !!          in the grid.
  subroutine get_B_elements(gauss_idx, eps, curr_weight, quadblock_B)
    use mod_equilibrium, only: rho_field
    use mod_make_subblock, only: subblock, reset_positions, reset_factors

    !> current index in the Gaussian grid
    integer, intent(in)          :: gauss_idx
    !> current value for the scale factor epsilon
    real(dp), intent(in)         :: eps
    !> current weight in the Gaussian quadrature
    real(dp), intent(in)         :: curr_weight
    !> the B-quadblock for a particular grid interval
    complex(dp), intent(inout)   :: quadblock_B(dim_quadblock, dim_quadblock)

    complex(dp), allocatable     :: factors(:)
    integer, allocatable         :: positions(:, :)

    real(dp)                     :: rho

    rho = rho_field % rho0(gauss_idx)

    ! Quadratic * Quadratic
    call reset_factors(factors, 5)
    call reset_positions(positions, 5)
    ! B(1,1)
    factors(1) = 1.0d0 / eps
    positions(1, :) = [1, 1]
    ! B(3,3)
    factors(2) = rho
    positions(2, :) = [3, 3]
    ! B(4,4)
    factors(3) = rho / eps
    positions(3, :) = [4, 4]
    ! B(5,5)
    factors(4) = rho / eps
    positions(4, :) = [5, 5]
    ! B(6,6)
    factors(5) = 1.0d0
    positions(5, :) = [6, 6]
    call subblock(quadblock_B, factors, positions, curr_weight, h_quadratic, h_quadratic)

    ! Cubic * Cubic
    call reset_factors(factors, 3)
    call reset_positions(positions, 3)
    ! B(2,2)
    factors(1) = rho / eps
    positions(1, :) = [2, 2]
    ! B(7,7)
    factors(2) = 1.0d0 / eps
    positions(2, :) = [7, 7]
    ! B(8,8)
    factors(3) = 1.0d0
    positions(3, :) = [8, 8]
    call subblock(quadblock_B, factors, positions, curr_weight, h_cubic, h_cubic)

    deallocate(factors)
    deallocate(positions)
  end subroutine get_B_elements



  !> Retrieves the A-matrix elements for a given quadblock corresponding
  !! to a particular point in the grid and a particular Gaussian weight.
  !! This routine is called <tt>n_gauss</tt> times for every grid interval.
  !! @note    When this subroutine gets called the basis functions defined at module scope
  !!          have already been recalculated and correspond to this particular value
  !!          in the grid.
  subroutine get_A_elements(gauss_idx, eps, d_eps_dr, curr_weight, quadblock_A)
    use mod_global_variables, only: ic, gamma_1
    use mod_equilibrium_params, only: k2, k3
    use mod_equilibrium, only: rho_field, T_field, B_field, v_field, grav_field, eta_field, &
                               rc_field, kappa_field
    use mod_make_subblock, only: subblock, reset_factors, reset_positions

    !> current index in the Gaussian grid
    integer, intent(in)       :: gauss_idx
    !> current value for the scale factor epsilon
    real(dp), intent(in)      :: eps
    !> current value for the derivative of epsilon
    real(dp), intent(in)      :: d_eps_dr
    !> current weight in the Gaussian quadrature
    real(dp), intent(in)      :: curr_weight
    !> the A-quadblock for a particular grid interval
    complex(dp), intent(inout):: quadblock_A(dim_quadblock, dim_quadblock)

    complex(dp), allocatable  :: factors(:)
    integer, allocatable      :: positions(:, :)

    real(dp)                  :: eps_inv
    real(dp)                  :: rho0, T0, B02, B03, B2_inv
    real(dp)                  :: v02, v03
    real(dp)                  :: tc_para, tc_perp, L0, eta, grav

    real(dp)                  :: drho0, drB02, dB02_r, dB03, dT0, dB02
    real(dp)                  :: drv02, dv02, dv03
    real(dp)                  :: dtc_perp_dT, dtc_perp_drho, dtc_perp_dB2
    real(dp)                  :: L_T, L_rho
    real(dp)                  :: deta_dT, deta, ddB03, ddB02

    eps_inv = 1.0d0 / eps

    !! Equilibrium quantities
    rho0 = rho_field % rho0(gauss_idx)
    drho0 = rho_field % d_rho0_dr(gauss_idx)

    T0 = T_field % T0(gauss_idx)
    dT0 = T_field % d_T0_dr(gauss_idx)

    B02 = B_field % B02(gauss_idx)
    dB02 = B_field % d_B02_dr(gauss_idx)
    B03 = B_field % B03(gauss_idx)
    dB03   = B_field % d_B03_dr(gauss_idx)
    B2_inv = 1.0d0 / (B_field % B0(gauss_idx))**2

    v02 = v_field % v02(gauss_idx)
    dv02 = v_field % d_v02_dr(gauss_idx)
    v03  = v_field % v03(gauss_idx)
    dv03 = v_field % d_v03_dr(gauss_idx)

    grav = grav_field % grav(gauss_idx)

    eta = eta_field % eta(gauss_idx)
    deta_dT = eta_field % d_eta_dT(gauss_idx)
    ddB02 = eta_field % dd_B02_dr(gauss_idx)
    ddB03 = eta_field % dd_B03_dr(gauss_idx)
    deta = eta_field % d_eta_dr(gauss_idx) + deta_dT * dT0

    L0 = rc_field % heat_loss(gauss_idx)
    L_T = rc_field % d_L_dT(gauss_idx)
    L_rho = rc_field % d_L_drho(gauss_idx)

    tc_para = kappa_field % kappa_para(gauss_idx)
    tc_perp = kappa_field % kappa_perp(gauss_idx)
    dtc_perp_drho = kappa_field % d_kappa_perp_drho(gauss_idx)
    dtc_perp_dT = kappa_field % d_kappa_perp_dT(gauss_idx)
    dtc_perp_dB2 = kappa_field % d_kappa_perp_dB2(gauss_idx)

    !! Calculate derivatives eps*B02, B02/eps, eps*V02
    drB02 = d_eps_dr * B02 + eps * dB02
    dB02_r = eps_inv * dB02 - d_eps_dr * eps_inv**2 * B02
    drv02 = d_eps_dr * v02 + eps * dv02


    !! Setup of matrix elements

    ! Quadratic * Quadratic
    call reset_factors(factors, 19)
    call reset_positions(positions, 19)
    ! A(1, 1)
    factors(1) = eps_inv * (eps_inv * v02 * k2 + v03 * k3)
    positions(1, :) = [1, 1]
    ! A(1, 3)
    factors(2) = eps_inv * rho0 * k2
    positions(2, :) = [1, 3]
    ! A(1, 4)
    factors(3) = eps_inv * rho0 * k3
    positions(3, :) = [1, 4]
    ! A(3, 1)
    factors(4) = eps_inv**2 * T0 * k2
    positions(4, :) = [3, 1]
    ! A(3, 3)
    factors(5) = rho0 * (eps_inv * k2 * v02 + k3 * v03)
    positions(5, :) = [3, 3]
    ! A(3, 5)
    factors(6) = eps_inv**2 * rho0 * k2
    positions(6, :) = [3, 5]
    ! A(3, 6)
    factors(7) = -B03 * (k3**2 + k2**2 * eps_inv**2)
    positions(7, :) = [3, 6]
    ! A(4, 1)
    factors(8) = eps_inv * T0 * k3
    positions(8, :) = [4, 1]
    ! A(4, 4)
    factors(9) = eps_inv * rho0 * (eps_inv * k2 * v02 + k3 * v03)
    positions(9, :) = [4, 4]
    ! A(4, 5)
    factors(10) = eps_inv * rho0 * k3
    positions(10, :) = [4, 5]
    ! A(4, 6)
    factors(11) = B02 * (k3**2 + k2**2 * eps_inv**2)
    positions(11, :) = [4, 6]
    ! A(5, 1)
    factors(12) = ic * gamma_1 * eps_inv * (d_eps_dr * eps_inv * dT0 * dtc_perp_drho - L0 - rho0*L_rho)
    positions(12, :) = [5, 1]
    ! A(5, 3)
    factors(13) = gamma_1 * eps_inv * rho0 * T0 * k2
    positions(13, :) = [5, 3]
    ! A(5, 4)
    factors(14) = gamma_1 * eps_inv * rho0 * T0 * k3
    positions(14, :) = [5, 4]
    ! A(5, 5)
    factors(15) = -eps_inv * ic * gamma_1 * ( tc_perp * (d_eps_dr * eps_inv)**2 &
                                              + tc_perp * (k2**2 * eps_inv**2 + k3**2) &
                                              + rho0 * L_T - d_eps_dr * eps_inv * dT0 * dtc_perp_dT ) &
                  + eps_inv * rho0 * (eps_inv * k2 * v02 + k3 * v03) &
                  + ic * gamma_1 * eps_inv * deta_dT * ( dB02**2 + dB03**2 + 2.0d0*d_eps_dr * eps_inv * B02 * dB02 &
                                                        + (d_eps_dr * eps_inv * B02)**2  )
    if (thermal_conduction) then
      factors(15) = factors(15) - eps_inv * ic * gamma_1 * (tc_para-tc_perp) * B2_inv * (k2 * eps_inv * B02 + k3 * B03)**2
    end if
    positions(15, :) = [5, 5]
    ! A(5, 6)
    factors(16) = 2.0d0 * ic * gamma_1 * ( (dT0 * d_eps_dr * eps_inv**2 * (eps * B02 * k3 - B03 * k2) * dtc_perp_dB2) &
                                          + eps_inv * k2 * eta * ddB03 &
                                          - k3 * eta * ddB02 &
                                          + 2 * k3 * (d_eps_dr * eps_inv)**2 * eta * B02 &
                                          + eps_inv * k2 * dB03 * deta &
                                          - k3 * dB02 * deta &
                                          - k3 * d_eps_dr * eps_inv * B02 * deta )
    positions(16, :) = [5, 6]
    ! A(6, 3)
    factors(17) = -B03
    positions(17, :) = [6, 3]
    ! A(6, 4)
    factors(18) = eps_inv * B02
    positions(18, :) = [6, 4]
    ! A(6, 6)
    factors(19) = eps_inv * k2 * v02 + k3 * v03 - ic * eta * (eps_inv**2 * k2**2 + k3**2)
    positions(19, :) = [6, 6]
    call subblock(quadblock_A, factors, positions, curr_weight, h_quadratic, h_quadratic)


    ! Quadratic * Cubic
    call reset_factors(factors, 10)
    call reset_positions(positions, 10)
    ! A(1, 2)
    factors(1) = -eps_inv * drho0
    positions(1, :) = [1, 2]
    ! A(3, 2)
    factors(2) = -eps_inv**2 * rho0 * drv02
    positions(2, :) = [3, 2]
    ! A(3, 7)
    factors(3) = eps_inv**2 * drB02 * k3
    positions(3, :) = [3, 7]
    ! A(3, 8)
    factors(4) = -eps_inv**2 * drB02 * k2
    positions(4, :) = [3, 8]
    ! A(4, 2)
    factors(5) = -eps_inv * rho0 * dv03
    positions(5, :) = [4, 2]
    ! A(4, 7)
    factors(6) = eps_inv * dB03 * k3
    positions(6, :) = [4, 7]
    ! A(4, 8)
    factors(7) = -eps_inv * dB03 * k2
    positions(7, :) = [4, 8]
    ! A(5, 2)
    factors(8) = -eps_inv * rho0 * dT0
    positions(8, :) = [5, 2]
    ! A(5, 7)
    factors(9) = ic * gamma_1 * eps_inv * (-2.0d0 * eta * k3**2 * dB03 - 2.0d0 * eta * k2 * k3 * eps_inv**2 * drB02)
    if (thermal_conduction) then
      factors(9) = factors(9) + ic * gamma_1 * eps_inv * ( (tc_para - tc_perp) * B2_inv * dT0 &
                                                            * (k2 * k3 * B02 * eps_inv + k3**2 * B03) )
    end if
    positions(9, :) = [5, 7]
    ! A(5, 8)
    factors(10) = -ic * gamma_1 * eps_inv * (-2.0d0 * eta * k2 * k3 * dB03 - 2.0d0 * eta * k2**2 * eps_inv**2 * drB02)
    if (thermal_conduction) then
      factors(10) = factors(10) - ic * gamma_1 * eps_inv * ( (tc_para - tc_perp) * B2_inv * dT0 &
                                                              * (k2**2 * B02 * eps_inv + k2 * k3 * B03) )
    end if
    positions(10, :) = [5, 8]
    call subblock(quadblock_A, factors, positions, curr_weight, h_quadratic, h_cubic)


    ! Quadratic * d(Cubic)/dr
    call reset_factors(factors, 10)
    call reset_positions(positions, 10)
    ! A(1, 2)
    factors(1) = -eps_inv * rho0
    positions(1, :) = [1, 2]
    ! A(3, 7)
    factors(2) = eps_inv**2 * B03 * k2
    positions(2, :) = [3, 7]
    ! A(3, 8)
    factors(3) = B03 * k3
    positions(3, :) = [3, 8]
    ! A(4, 7)
    factors(4) = -eps_inv**2 * B02 * k2
    positions(4, :) = [4, 7]
    ! A(4, 8)
    factors(5) = -B02 * k3
    positions(5, :) = [4, 8]
    ! A(5, 2)
    factors(6) = -gamma_1 * eps_inv * rho0 * T0
    positions(6, :) = [5, 2]
    ! A(5, 7)
    factors(7) = 2.0d0 * ic * gamma_1 * eps_inv * ( dT0 * B03 * d_eps_dr * eps_inv * dtc_perp_dB2 &
                                                    - eta * ddB03 - deta * dB03 )
    positions(7, :) = [5, 7]
    ! A(5, 8)
    factors(8) = -2.0d0 * ic * gamma_1 * ( dT0 * d_eps_dr * eps_inv * B02 * dtc_perp_dB2 &
                                          - eta * ddB02 &
                                          + 2.0d0 * (d_eps_dr * eps_inv)**2 * eta * B02 &
                                          - dB02 * deta &
                                          - d_eps_dr * eps_inv * B02 * deta )
    positions(8, :) = [5, 8]
    ! A(6, 7)
    factors(9) = -eps_inv * v02 + ic * eta * eps_inv**2 * k2
    positions(9, :) = [6, 7]
    ! A(6, 8)
    factors(10) = -v03 + ic * eta * k3
    positions(10, :) = [6, 8]
    call subblock(quadblock_A, factors, positions, curr_weight, h_quadratic, dh_cubic_dr)


    ! Cubic * Quadratic
    call reset_factors(factors, 7)
    call reset_positions(positions, 7)
    ! A(2, 1)
    factors(1) = -d_eps_dr * eps_inv**2 * v02**2 + eps_inv * grav
    positions(1, :) = [2, 1]
    ! A(2, 3)
    factors(2) = -2.0d0 * d_eps_dr * eps_inv * rho0 * v02
    positions(2, :) = [2, 3]
    ! A(2, 6)
    factors(3) = eps_inv * drB02 * k3 - eps * k3 * dB02_r
    positions(3, :) = [2, 6]
    ! A(7, 5)
    factors(4) = ic * eps_inv * deta_dT * dB03
    positions(4, :) = [7, 5]
    ! A(7, 6)
    factors(5) = ic * deta * eps_inv * k2
    positions(5, :) = [7, 6]
    ! A(8, 5)
    factors(6) = -ic * eps_inv**2 * deta_dT * drB02
    positions(6, :) = [8, 5]
    ! A(8, 6)
    factors(7) = -ic * ( eta * d_eps_dr * eps_inv - deta) * k3
    positions(7, :) = [8, 6]
    call subblock(quadblock_A, factors, positions, curr_weight, h_cubic, h_quadratic)


    ! d(Cubic)/dr * Quadratic
    call reset_factors(factors, 5)
    call reset_positions(positions, 5)
    ! A(2, 1)
    factors(1) = -T0 * eps_inv
    positions(1, :) = [2, 1]
    ! A(2, 5)
    factors(2) = -rho0 * eps_inv
    positions(2, :) = [2, 5]
    ! A(2, 6)
    factors(3) = eps_inv * k2 * B03 - k3 * B02
    positions(3, :) = [2, 6]
    ! A(7, 6)
    factors(4) = ic * eta * eps_inv * k2
    positions(4, :) = [7, 6]
    ! A(8, 6)
    factors(5) = ic * eta * k3
    positions(5, :) = [8, 6]
    call subblock(quadblock_A, factors, positions, curr_weight, dh_cubic_dr, h_quadratic)


    ! Cubic * Cubic
    call reset_factors(factors, 9)
    call reset_positions(positions, 9)
    ! A(2, 2)
    factors(1) = eps_inv * (eps_inv * v02 * k2 + v03 * k3) * rho0
    positions(1, :) = [2, 2]
    ! A(2, 7)
    factors(2) = -B03 * eps_inv * k3**2 - B02 * eps_inv**2 * k2 * k3
    positions(2, :) = [2, 7]
    ! A(2, 8)
    factors(3) = B03 * eps_inv * k2 * k3 + B02 * eps_inv**2 * k2**2
    positions(3, :) = [2, 8]
    ! A(7, 2)
    factors(4) = -eps_inv * B03
    positions(4, :) = [7, 2]
    ! A(7, 7)
    factors(5) = eps_inv * k3 * (v03 - ic * eta * k3)
    positions(5, :) = [7, 7]
    ! A(7, 8)
    factors(6) = -eps_inv * k2 * (v03 - ic * eta * k3)
    positions(6, :) = [7, 8]
    ! A(8, 2)
    factors(7) = eps_inv * B02
    positions(7, :) = [8, 2]
    ! A(8, 7)
    factors(8) = -eps_inv * k3 * (v02 - ic * eta * eps_inv * k2)
    positions(8, :) = [8, 7]
    ! A(8, 8)
    factors(9) = eps_inv * k2 * (v02 - ic * eta * eps_inv * k2)
    positions(9, :) = [8, 8]
    call subblock(quadblock_A, factors, positions, curr_weight, h_cubic, h_cubic)


    ! Cubic * d(Cubic)/dr
    call reset_factors(factors, 3)
    call reset_positions(positions, 3)
    ! A(2, 8)
    factors(1) = -2.0d0 * B02 * d_eps_dr * eps_inv
    positions(1, :) = [2, 8]
    ! A(7, 7)
    factors(2) = -ic * deta * eps_inv
    positions(2, :) = [7, 7]
    ! A(8, 8)
    factors(3) = ic * (eta * d_eps_dr * eps_inv - deta)
    positions(3, :) = [8, 8]
    call subblock(quadblock_A, factors, positions, curr_weight, h_cubic, dh_cubic_dr)


    ! d(Cubic)/dr * d(Cubic)/dr
    call reset_factors(factors, 4)
    call reset_positions(positions, 4)
    ! A(2, 7)
    factors(1) = -B03 * eps_inv
    positions(1, :) = [2, 7]
    ! A(2, 8)
    factors(2) = B02
    positions(2, :) = [2, 8]
    ! A(7, 7)
    factors(3) = -ic * eta * eps_inv
    positions(3, :) = [7, 7]
    ! A(8, 8)
    factors(4) = -ic * eta
    positions(4, :) = [8, 8]
    call subblock(quadblock_A, factors, positions, curr_weight, dh_cubic_dr, dh_cubic_dr)


    ! d(Quadratic)/dr * Quadratic
    call reset_factors(factors, 3)
    call reset_positions(positions, 3)
    ! A(5, 1)
    factors(1) = -ic * gamma_1 * eps_inv * dT0 * dtc_perp_drho
    positions(1, :) = [5, 1]
    ! A(5, 5)
    factors(2) = ic * gamma_1 * eps_inv * (d_eps_dr * eps_inv * tc_perp - dT0 * dtc_perp_dT)
    positions(2, :) = [5, 5]
    ! A(5, 6)
    factors(3) = -2.0d0 * ic * gamma_1 * eps_inv * ( dT0 * (eps * B02 * k3 - B03 * k2) * dtc_perp_dB2 &
                                                    - eta * k2 * dB03 + eta * k3 * drB02 )
    positions(3, :) = [5, 6]
    call subblock(quadblock_A, factors, positions, curr_weight, dh_quadratic_dr, h_quadratic)


    ! Quadratic * d(Quadratic)/dr
    call reset_factors(factors, 1)
    call reset_positions(positions, 1)
    ! A(5, 5)
    factors(1) = ic * gamma_1 * d_eps_dr * eps_inv**2 * tc_perp
    positions(1, :) = [5, 5]
    call subblock(quadblock_A, factors, positions, curr_weight, h_quadratic, dh_quadratic_dr)


    ! d(Quadratic)/dr * d(Quadratic)/dr
    call reset_factors(factors, 1)
    call reset_positions(positions, 1)
    ! A(5, 5)
    factors(1) = -ic * gamma_1 * eps_inv * tc_perp
    positions(1, :) = [5, 5]
    call subblock(quadblock_A, factors, positions, curr_weight, dh_quadratic_dr, dh_quadratic_dr)


    ! d(Quadratic)/dr * d(Cubic)/dr
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(5, 7)
    factors(1) = -2.0d0 * ic * gamma_1 * eps_inv * (dT0 * B03 * dtc_perp_dB2 + eta * dB03)
    positions(1, :) = [5, 7]
    ! A(5, 8)
    factors(2) = 2.0d0 * ic * gamma_1 * (dT0 * B02 * dtc_perp_dB2 + eta * eps_inv * drB02)
    positions(2, :) = [5, 8]
    call subblock(quadblock_A, factors, positions, curr_weight, dh_quadratic_dr, dh_cubic_dr)

    deallocate(factors)
    deallocate(positions)
  end subroutine get_A_elements

end module mod_matrix_creation
