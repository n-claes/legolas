module mod_radiative_cooling
  use mod_global_variables
  implicit none

  private

  !> Radiative cooling table containing temperatures
  real(dp), allocatable         :: interp_table_T(:)
  !> Radiative cooling table containing luminosities
  real(dp), allocatable         :: interp_table_L(:)
  !> Radiative cooling table containing derivative luminosity to temperature
  real(dp), allocatable         :: interp_table_dLdT(:)
  !> Log10 of minimum temperature in the interpolated cooling curve.
  real(dp)                      :: lgmin_T
  !> Log10 of maximum temperature in the interpolated cooling curve.
  real(dp)                      :: lgmax_T
  !> Log10 of the stepsize in the interpolated cooling curve.
  real(dp)                      :: lgstep
  !> Maximum temperature of interpolated cooling curve.
  real(dp)                      :: max_T
  !> Minimum temperature of interpolated cooling curve.
  real(dp)                      :: min_T

  public  :: initialise_radiative_cooling
  public  :: get_Lambda
  public  :: get_dLambdadT
  public  :: radiative_cooling_clean

contains

  !> Initialises the radiative cooling module, depending on the cooling
  !  curve specified.
  subroutine initialise_radiative_cooling()
    use mod_cooling_curves

    real(dp), allocatable       :: table_T(:), table_L(:)
    integer                     :: ntable

    allocate(interp_table_T(ncool))
    allocate(interp_table_L(ncool))
    allocate(interp_table_dLdT(ncool))

    ! Initialise interpolated tables to zero
    interp_table_T    = 0.0d0
    interp_table_L    = 0.0d0
    interp_table_dLdT = 0.0d0

    ! Read cooling curve
    select case(cooling_curve)

    case('JCcorona')
       write(*,*) "Use Colgan & Feldman (2008) cooling curve"

       ntable = n_JCcorona

       allocate(table_T(1:ntable))
       allocate(table_L(1:ntable))
       table_T(1:ntable) = t_JCcorona(1:n_JCcorona)
       table_L(1:ntable) = l_JCcorona(1:n_JCcorona)

    case('DM')
       write(*,*) "Use Delgano & McCray (1972) cooling curve"

       ntable = n_DM

       allocate(table_T(1:ntable))
       allocate(table_L(1:ntable))
       table_T(1:ntable) = t_DM(1:n_DM)
       table_L(1:ntable) = l_DM(1:n_DM)

    case('MLsolar')
       write(*,*) "Use Mellema & Lundqvist (2002) cooling curve ", &
                  "for solar metallicity"

       ntable = n_MLsolar

       allocate(table_T(1:ntable))
       allocate(table_L(1:ntable))
       table_T(1:ntable) = t_MLsolar(1:n_MLsolar)
       table_L(1:ntable) = l_MLsolar(1:n_MLsolar)

    case('SPEX')
       write(*,*) "Use SPEX cooling curve (Schure et al. 2009) ", &
                  "for solar metallicity"

       ntable = n_SPEX

       allocate(table_T(1:ntable))
       allocate(table_L(1:ntable))
       table_T(1:ntable) = t_SPEX(1:n_SPEX)
       table_L(1:ntable) = l_SPEX(1:n_SPEX) + log10(nenh_SPEX(1:n_SPEX))

    case('SPEX_DM')
          write(*,*) "Use SPEX cooling curve for solar metallicity ", &
                     "above 10^4 K, Schure et al. (2009)."
          write(*,*) "At lower temperatures, use Dalgarno & McCray (1972), ", &
                     "with a pre-set ionization fraction of 10^-3."

       ntable = n_SPEX + n_DM_2 - 6

       allocate(table_T(1:ntable))
       allocate(table_L(1:ntable))
       table_T(1:n_DM_2-1) = t_DM_2(1:n_DM_2-1)
       table_L(1:n_DM_2-1) = L_DM_2(1:n_DM_2-1)
       table_T(n_DM_2:ntable) = t_SPEX(6:n_SPEX)
       table_L(n_DM_2:ntable) = l_SPEX(6:n_SPEX) + log10(nenh_SPEX(6:n_SPEX))
    case default
      write(*, *) "Cooling curve not defined correctly."
      write(*, *) "Currently set on:   ", cooling_curve
      stop
    end select

    ! Interpolate cooling curves
    call interpolate_cooling_curve(ntable, table_T, table_L)

    deallocate(table_T)
    deallocate(table_L)

  end subroutine initialise_radiative_cooling


  !> Interpolates the luminosity based on the given temperatures T0.
  !! @param T0      Equilibrium temperatures, in K.
  !!                (array, length = 4*gridpts)
  !! @param lambda  Interpolated lambda for every T0
  !!                (array, length = 4*gridpts)
  subroutine get_Lambda(T0, lambda)
    real(dp), intent(in)  :: T0(4*gridpts)
    real(dp), intent(out) :: lambda(4*gridpts)

    integer               :: idx, i

    do i = 1, 4*gridpts
      idx = int( (log10(T0(i)) - lgmin_T) / lgstep ) + 1
      lambda(i) = interp_table_L(idx) + (T0(i) - interp_table_T(idx)) &
                  * (interp_table_L(idx + 1) - interp_table_L(idx))   &
                  / (interp_table_T(idx + 1) - interp_table_T(idx))
    end do

  end subroutine get_Lambda


  !> Interpolates the derivative of the cooling curve in T0.
  !! @param T0    Equilibrium temperatures, in K.
  !!              (array, length=4*gridpts)
  !! @param dLambdadT  Derivative of cooling curve with respect to temperature,
  !!                   evaluated for every T0.
  !!                   (array, length=4*gridpts)
  subroutine get_dLambdadT(T0, dLambdadT)
    real(dp), intent(in)  :: T0(4*gridpts)
    real(dp), intent(out) :: dLambdadT(4*gridpts)
    integer               :: idx, i

    do i = 1, 4*gridpts
      idx     = int( (log10(T0(i)) - lgmin_T) / lgstep ) + 1
      dLambdadT(i) = interp_table_dLdT(idx) + (T0(i) - interp_table_T(idx)) &
                     * (interp_table_dLdT(idx+1) - interp_table_dLdT(idx))  &
                     / (interp_table_T(idx+1)    - interp_table_T(idx))
    end do

  end subroutine get_dLambdadT


  !> Interpolates the rudimentary cooling curves using ncool points.
  !! A second-order polynomial interpolation is used except near sharp jumps.
  !! @param ntable    number of entries in the cooling table.
  !! @param table_T   Temperature entries of cooling table.
  !! @param table_L   Luminosity entries of cooling table.
  subroutine interpolate_cooling_curve(ntable, table_T, table_L)
    use mod_physical_constants

    real(dp), intent(in)  :: table_T(:), table_L(:)
    integer, intent(in)   :: ntable
    real(dp)              :: fact1, fact2, fact3, dL1, dL2, ratt
    integer               :: i, j
    logical               :: jump

    max_T = table_T(ntable)
    min_T = table_T(1)
    ratt     = (max_T-min_T) / dble(ncool-1)

    interp_table_T(1) = min_T
    interp_table_L(1) = table_L(1)

    interp_table_T(ncool) = max_T
    interp_table_L(ncool) = table_L(ntable)

    do i=2, ncool        ! loop to create one table
      interp_table_T(i) = interp_table_T(i-1) + ratt
      do j=1,ntable-1    ! loop to create one spot on a table
      ! Second order polynomial interpolation, except at the outer edge,
      ! or in case of a large jump.
        if (interp_table_T(i) < table_T(j+1)) then
           if (j.eq. ntable-1 ) then
             fact1 = (interp_table_T(i)-table_T(j+1))      &
                   /(table_T(j)-table_T(j+1))

             fact2 = (interp_table_T(i)-table_T(j))        &
                   /(table_T(j+1)-table_T(j))

             interp_table_L(i) = table_L(j)*fact1 + table_L(j+1)*fact2
             exit
           else
             dL1 = table_L(j+1)-table_L(j)
             dL2 = table_L(j+2)-table_L(j+1)
             jump =(max(dabs(dL1),dabs(dL2)) > 2*min(dabs(dL1),dabs(dL2)))
           endif

           if( jump ) then
             fact1 = (interp_table_T(i)-table_T(j+1))      &
                   /(table_T(j)-table_T(j+1))

             fact2 = (interp_table_T(i)-table_T(j))        &
                   /(table_T(j+1)-table_T(j))

             interp_table_L(i) = table_L(j)*fact1 + table_L(j+1)*fact2
             exit
           else
             fact1 = ((interp_table_T(i)-table_T(j+1))     &
                   * (interp_table_T(i)-table_T(j+2)))     &
                   / ((table_T(j)-table_T(j+1))            &
                   * (table_T(j)-table_T(j+2)))

             fact2 = ((interp_table_T(i)-table_T(j))       &
                   * (interp_table_T(i)-table_T(j+2)))     &
                   / ((table_T(j+1)-table_T(j))            &
                   * (table_T(j+1)-table_T(j+2)))

             fact3 = ((interp_table_T(i)-table_T(j))       &
                   * (interp_table_T(i)-table_T(j+1)))     &
                   / ((table_T(j+2)-table_T(j))            &
                   * (table_T(j+2)-table_T(j+1)))

             interp_table_L(i) = table_L(j)*fact1 + table_L(j+1)*fact2 &
                      + table_L(j+2)*fact3
             exit
           endif
        endif
      enddo  ! end loop to find create one spot on a table
    enddo    ! end loop to create one table

    ! Go from logarithmic to actual values.
    interp_table_T(1:ncool) = 10.0d0**interp_table_T(1:ncool)
    interp_table_L(1:ncool) = 10.0d0**interp_table_L(1:ncool)

    ! Normalise values
    interp_table_T(1:ncool) = interp_table_T(1:ncool) / unit_temperature
    interp_table_L(1:ncool) = interp_table_L(1:ncool) / unit_luminosity

    lgmin_T = dlog10(min_T)
    lgmax_T = dlog10(max_T)
    lgstep = (lgmax_T-lgmin_T) * 1.d0 / (ncool-1)

    interp_table_dLdT(1)     = (interp_table_L(2) - interp_table_L(1)) &
                             / (interp_table_T(2) - interp_table_T(1))
    interp_table_dLdT(ncool) = (interp_table_L(ncool) - interp_table_L(ncool-1)) &
                             / (interp_table_T(ncool) - interp_table_T(ncool-1))

    do i=2,ncool-1
      interp_table_dLdT(i) = (interp_table_L(i+1) - interp_table_L(i-1)) &
                           / (interp_table_T(i+1) - interp_table_T(i-1))
    enddo

  end subroutine interpolate_cooling_curve

  !> Deallocates arrays defined in the radiative cooling module.
  subroutine radiative_cooling_clean()
    deallocate(interp_table_T)
    deallocate(interp_table_L)
    deallocate(interp_table_dLdT)
  end subroutine radiative_cooling_clean





end module mod_radiative_cooling
