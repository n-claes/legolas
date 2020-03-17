!
! MODULE: mod_radiative_cooling
!
!> @author
!> Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Module to calculate the radiative cooling contributions.
!
module mod_radiative_cooling
  use mod_global_variables, only: dp, ncool
  implicit none

  private

  !> Radiative cooling table containing temperatures
  real(dp), allocatable :: interp_table_T(:)
  !> Radiative cooling table containing luminosities
  real(dp), allocatable :: interp_table_L(:)
  !> Radiative cooling table containing derivative luminosity to temperature
  real(dp), allocatable :: interp_table_dLdT(:)
  !> Log10 of minimum temperature in the interpolated cooling curve.
  real(dp)              :: lgmin_T
  !> Log10 of maximum temperature in the interpolated cooling curve.
  real(dp)              :: lgmax_T
  !> Log10 of the stepsize in the interpolated cooling curve.
  real(dp)              :: lgstep
  !> Maximum temperature of interpolated cooling curve.
  real(dp)              :: max_T
  !> Minimum temperature of interpolated cooling curve.
  real(dp)              :: min_T
  !> Use an interpolated cooling curve or not
  logical, save         :: interpolated_curve = .true.

  public  :: initialise_radiative_cooling
  public  :: set_radiative_cooling_values
  public  :: radiative_cooling_clean

contains

  !> Initialises the radiative cooling module, depending on the cooling
  !! curve specified. Interpolates the data from mod_cooling_curves using
  !! ncool gridpoints.
  subroutine initialise_radiative_cooling()
    use mod_global_variables, only: cooling_curve
    use mod_cooling_curves

    real(dp), allocatable :: table_T(:), table_L(:)
    integer               :: ntable

    select case(cooling_curve)
    case('jc_corona')
       ntable = n_jc_corona
       allocate(table_T(ntable))
       allocate(table_L(ntable))
       table_T = t_jc_corona
       table_L = l_jc_corona

    case('dalgarno')
       ntable = n_dalgarno
       allocate(table_T(ntable))
       allocate(table_L(ntable))
       table_T = t_dalgarno
       table_L = l_dalgarno

    case('ml_solar')
       ntable = n_ml_solar
       allocate(table_T(ntable))
       allocate(table_L(ntable))
       table_T = t_ml_solar
       table_L = l_ml_solar

    case('spex')
       ntable = n_spex
       allocate(table_T(ntable))
       allocate(table_L(ntable))
       table_T = t_spex
       table_l = l_spex

    case('spex_dalgarno')
       ntable = n_spex + n_dalgarno2 - 6
       allocate(table_T(ntable))
       allocate(table_L(ntable))
       table_T(1:n_dalgarno2-1) = t_dalgarno2(1:n_dalgarno2-1)
       table_L(1:n_dalgarno2-1) = l_dalgarno2(1:n_dalgarno2-1)
       table_T(n_dalgarno2:ntable) = t_spex(6:n_spex)
       table_L(n_dalgarno2:ntable) = l_spex(6:n_SPEX) + log10(n_spex_enh(6:n_spex))

    case('rosner')
      interpolated_curve = .false.

    case default
      write(*, *) "Unknown cooling curve '", cooling_curve, "' provided!"
      error stop
    end select

    if (interpolated_curve) then
      allocate(interp_table_T(ncool))
      allocate(interp_table_L(ncool))
      allocate(interp_table_dLdT(ncool))

      ! Initialise interpolated tables to zero
      interp_table_T    = 0.0d0
      interp_table_L    = 0.0d0
      interp_table_dLdT = 0.0d0

      call interpolate_cooling_curve(ntable, table_T, table_L)

      deallocate(table_T)
      deallocate(table_L)
    end if

  end subroutine initialise_radiative_cooling


  subroutine set_radiative_cooling_values(rho_field, T_field, rc_field)
    use mod_types, only: density_type, temperature_type, cooling_type
    use mod_global_variables, only: gauss_gridpts, cooling_curve
    use mod_cooling_curves, only: get_rosner_cooling

    type(density_type), intent(in)      :: rho_field
    type(temperature_type), intent(in)  :: T_field
    type(cooling_type), intent(inout)   :: rc_field

    real(dp)    :: lambda_T(gauss_gridpts)
    real(dp)    :: d_lambda_dT(gauss_gridpts)
    real(dp)    :: T0(gauss_gridpts)
    integer     :: idx, i

    T0 = T_field % T0

    if (interpolated_curve) then
      do i = 1, gauss_gridpts
          if (T0(i) <= min_T) then
            ! No cooling if T0 below lower limit of cooling curve
            lambda_T(i) = 0.0d0
            d_lambda_dT(i) = 0.0d0
          else if (T0(i) >= max_T) then
            ! Assume Bremmstrahlung sqrt(T/Tmax) if T is above upper limit of cooling curve
            lambda_T(i) = interp_table_L(ncool) * sqrt(T0(i) / max_T)
            d_lambda_dT(i) = 0.5d0 * interp_table_L(ncool) / sqrt(T0(i) * max_T)
          else
            ! Interpolate values from the cooling curves
            idx = int( (log10(T0(i)) - lgmin_T) / lgstep ) + 1
            lambda_T(i) = interp_table_L(idx) + (T0(i) - interp_table_T(idx)) &
                          * (interp_table_L(idx + 1) - interp_table_L(idx))   &
                          / (interp_table_T(idx + 1) - interp_table_T(idx))
            d_lambda_dT(i) = interp_table_dLdT(idx) + (T0(i) - interp_table_T(idx)) &
                             * (interp_table_dLdT(idx+1) - interp_table_dLdT(idx))  &
                             / (interp_table_T(idx+1) - interp_table_T(idx))
          end if
      end do
    else
      ! In this case an analytical cooling curve is used
      select case(cooling_curve)
      case('rosner')
        call get_rosner_cooling(T0, lambda_T, d_lambda_dT)

      case default
        write(*, *) "Unknown cooling curve '", cooling_curve, "' provided!"
        error stop
      end select
    end if

    ! dL/dT = rho0 * d_lambda_dT where lambda_T equals the cooling curve
    rc_field % d_L_dT = (rho_field % rho0) * d_lambda_dT
    ! dL/drho = lambda(T)
    rc_field % d_L_drho = lambda_T

  end subroutine set_radiative_cooling_values


  !> Interpolates the rudimentary cooling curves using ncool points.
  !! A second-order polynomial interpolation is used except near sharp jumps.
  !! @param[in] ntable    number of entries in the cooling table.
  !! @param[in] table_T   Temperature entries of cooling table.
  !! @param[in] table_L   Luminosity entries of cooling table.
  subroutine interpolate_cooling_curve(ntable, table_T, table_L)
    use mod_units, only: unit_temperature, unit_luminosity

    integer, intent(in)   :: ntable
    real(dp), intent(in)  :: table_T(:), table_L(:)
    real(dp)              :: fact1, fact2, fact3, dL1, dL2, ratt
    integer               :: i, j
    logical               :: jump

    max_T = table_T(ntable)
    min_T = table_T(1)
    ratt  = (max_T-min_T) / dble(ncool-1)

    interp_table_T(1) = min_T
    interp_table_L(1) = table_L(1)

    interp_table_T(ncool) = max_T
    interp_table_L(ncool) = table_L(ntable)

    do i=2, ncool        ! loop to create one table
      interp_table_T(i) = interp_table_T(i-1) + ratt
      do j=1, ntable-1    ! loop to create one spot on a table
      ! Second order polynomial interpolation, except at the outer edge,
      ! or in case of a large jump.
        if (interp_table_T(i) < table_T(j+1)) then
           if (j == ntable-1 ) then
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

           if (jump) then
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
      enddo  ! end loop to create one spot on a table
    enddo    ! end loop to create one table

    ! Go from logarithmic to actual values.
    interp_table_T(1:ncool) = 10.0d0**interp_table_T(1:ncool)
    interp_table_L(1:ncool) = 10.0d0**interp_table_L(1:ncool)

    ! Normalise values
    interp_table_T(1:ncool) = interp_table_T(1:ncool) / unit_temperature
    interp_table_L(1:ncool) = interp_table_L(1:ncool) / unit_luminosity

    min_T = interp_table_T(1)
    max_T = interp_table_T(ncool)

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
    if (interpolated_curve) then
      deallocate(interp_table_T)
      deallocate(interp_table_L)
      deallocate(interp_table_dLdT)
    end if
  end subroutine radiative_cooling_clean

end module mod_radiative_cooling
