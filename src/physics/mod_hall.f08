! =============================================================================
!> Module containing Hall-related routines.
!! Sets the Hall and electron inertia factors based on normalisations and specified profiles.
module mod_hall
  use mod_global_variables, only: dp
  use mod_settings, only: settings_t

implicit none

private

public  :: set_hall_factors

contains

  !> Retrieves the normalised Hall factor as described by Porth et al. (2014),
  !! with a dropoff at the boundary, if desired. Additionally, defines the
  !! electron inertia factor if included, with a dropoff profile, if desired.
  subroutine set_hall_factors(settings, hall_field)
    use mod_grid, only: grid_gauss
    use mod_physical_constants, only: dpi, mp_cgs, ec_cgs, me_cgs
    use mod_types, only: hall_type

    type(settings_t), intent(in) :: settings
    type(hall_type), intent(inout)  :: hall_field

    real(dp) :: unit_velocity, unit_length, unit_magneticfield
    real(dp)  :: sleft, sright, width, hallval, inertiaval, edge_dist
    real(dp)  :: x, shift, stretch, shift2, stretch2
    integer   :: i, gauss_gridpts

    gauss_gridpts = settings%grid%get_gauss_gridpts()
    width = settings%physics%dropoff_width
    edge_dist = settings%physics%dropoff_edge_dist

    unit_velocity = settings%units%get_unit_velocity()
    unit_length = settings%units%get_unit_length()
    unit_magneticfield = settings%units%get_unit_magneticfield()
    hallval = (mp_cgs * unit_velocity) / (ec_cgs * unit_length * unit_magneticfield)
    inertiaval = ( &
      mp_cgs * me_cgs * unit_velocity**2 &
      / (ec_cgs * unit_length * unit_magneticfield)**2 &
    )

    sleft = grid_gauss(1) + 0.5d0 * width + edge_dist
    sright = grid_gauss(gauss_gridpts) - edge_dist - 0.5d0 * width

    if (settings%physics%hall%use_dropoff) then
      shift = hallval * tanh(-dpi) / (tanh(-dpi) - tanh(dpi))
      stretch = hallval / (tanh(dpi) - tanh(-dpi))

      do i = 1, gauss_gridpts
        x = grid_gauss(i)
        if (sleft - 0.5d0*width <= x .and. x <= sleft + 0.5d0*width) then
          hall_field % hallfactor(i) = shift + stretch * tanh(2.0d0 * dpi * (x - sleft) / width)
        else if (sleft + 0.5d0*width < x .and. x < sright - 0.5d0*width) then
          hall_field % hallfactor(i) = hallval
        else if (sright - 0.5d0*width <= x .and. x <= sright + 0.5d0*width) then
          hall_field % hallfactor(i) = shift + stretch * tanh(2.0d0 * dpi * (sright - x) / width)
        else
          hall_field % hallfactor(i) = 0.0d0
        end if
      end do
    else
      hall_field % hallfactor = hallval
    end if

    if (settings%physics%hall%has_electron_inertia()) then
      if (settings%physics%hall%use_inertia_dropoff) then
        shift2 = inertiaval * tanh(-dpi) / (tanh(-dpi) - tanh(dpi))
        stretch2 = inertiaval / (tanh(dpi) - tanh(-dpi))

        do i = 1, gauss_gridpts
          x = grid_gauss(i)
          if (sleft - 0.5d0*width <= x .and. x <= sleft + 0.5d0*width) then
            hall_field % inertiafactor(i) = shift2 + stretch2 * tanh(2.0d0 * dpi * (x - sleft) / width)
          else if (sleft + 0.5d0*width < x .and. x < sright - 0.5d0*width) then
            hall_field % inertiafactor(i) = inertiaval
          else if (sright - 0.5d0*width <= x .and. x <= sright + 0.5d0*width) then
            hall_field % inertiafactor(i) = shift2 + stretch2 * tanh(2.0d0 * dpi * (sright - x) / width)
          else
            hall_field % inertiafactor(i) = 0.0d0
          end if
        end do
      else
        hall_field % inertiafactor = inertiaval
      end if
    end if
  end subroutine set_hall_factors

end module mod_hall
