module mod_setup_equilibrium
  implicit none

  real, allocatable         ::  rho_0(:)
  real, allocatable         ::  v_0(:, :, :)
  real, allocatable         ::  T_0(:)
  real, allocatable         ::  B_0(:, :, :)


contains

  subroutine initialise_equilibrium(grid, grid_gauss)
    use mod_global_variables

    double precision, intent(in)  :: grid(gridpts)
    double precision, intent(out) :: grid_gauss(4*gridpts)

    allocate(rho_0(gridpts))
    allocate(v_0(gridpts, gridpts, gridpts))
    allocate(T_0(gridpts))
    allocate(B_0(gridpts, gridpts, gridpts))

    call set_equilibrium(grid, grid_gauss)

  end subroutine initialise_equilibrium

  !> Sets the equilibrium on the nodes of the Gaussian quadrature
  subroutine set_equilibrium(grid, grid_gauss)
    use mod_global_variables

    double precision, intent(in)  :: grid(gridpts)
    double precision, intent(out) :: grid_gauss(4*gridpts)
    double precision              :: x_lo, x_hi, dx, xi(n_gauss)
    integer                       :: i, j, idx

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
    do i = 1, gridpts
      rho_0(i)     = 1.0d0
      v_0(i, i, i) = 1.0d0
      T_0(i)       = 1.0d0
      B_0(i, i, i) = 1.0d0
    end do

  end subroutine set_equilibrium

  subroutine grid_LEDA(SGRID)
    use mod_global_variables
    double precision, intent(in) :: SGRID(gridpts)
    double precision             :: SGI(4*gridpts)
    double precision :: SL, SU, ZDIF, ZA, ZB, ZC, ZS(4)
    integer :: I, J, NGINT, NI

    NGINT = gridpts - 1

    DO 90 NI=1,NGINT
       SL    = SGRID(NI)
       SU    = SGRID(NI+1)

! ... BERECHNUNG DER STUETZSTELLEN ...

       ZDIF  = SU - SL
       ZA    = .5 * ( SU + SL )
       ZB    = .43056815579702629 * ZDIF
       ZC    = .16999052179242813 * ZDIF
       ZS(1) = ZA - ZB
       ZS(2) = ZA - ZC
       ZS(3) = ZA + ZC
       ZS(4) = ZA + ZB

       DO 80 I=1,4
          J = (NI-1)*4+I
          SGI(J) = ZS(I)
 80    CONTINUE
    90 CONTINUE
    print*,SGI

  end subroutine grid_LEDA

  subroutine equilibrium_clean()
    deallocate(rho_0)
    deallocate(v_0)
    deallocate(T_0)
    deallocate(B_0)
  end subroutine equilibrium_clean


end module mod_setup_equilibrium
