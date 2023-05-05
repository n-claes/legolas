module mod_physics_utils
  use mod_global_variables, only: dp
  use mod_physical_constants, only: dpi
  use mod_settings, only: settings_t
  implicit none

  private

  public :: get_dropoff
  public :: get_dropoff_dr

contains


  real(dp) function get_dropoff(x, cte_value, settings)
    real(dp), intent(in) :: x
    real(dp), intent(in) :: cte_value
    type(settings_t), intent(in) :: settings
    real(dp) :: edge_dist, width, sleft, sright, shift, stretch

    edge_dist = settings%physics%dropoff_edge_dist
    width = settings%physics%dropoff_width
    sleft = settings%grid%get_grid_start() + 0.5_dp * width + edge_dist
    sright = settings%grid%get_grid_end() - 0.5_dp * width - edge_dist
    shift = cte_value * tanh(-dpi) / (tanh(-dpi) - tanh(dpi))
    stretch = cte_value / (tanh(dpi) - tanh(-dpi))

    if (sleft - 0.5_dp * width <= x .and. x <= sleft + 0.5_dp * width) then
      get_dropoff = shift + stretch * tanh(2.0_dp * dpi * (x - sleft) / width)
    else if (sleft + 0.5_dp * width < x .and. x < sright - 0.5_dp * width) then
      get_dropoff = cte_value
    else if (sright - 0.5_dp * width <= x .and. x <= sright + 0.5_dp * width) then
      get_dropoff = shift + stretch * tanh(2.0_dp * dpi * (sright - x) / width)
    else
      get_dropoff = 0.0_dp
    end if
  end function get_dropoff


  real(dp) function get_dropoff_dr(x, cte_value, settings)
    real(dp), intent(in) :: x
    real(dp), intent(in) :: cte_value
    type(settings_t), intent(in) :: settings
    real(dp) :: edge_dist, width, sleft, sright, shift, stretch

    edge_dist = settings%physics%dropoff_edge_dist
    width = settings%physics%dropoff_width
    sleft = settings%grid%get_grid_start() + 0.5_dp * width + edge_dist
    sright = settings%grid%get_grid_end() - 0.5_dp * width - edge_dist
    shift = cte_value * tanh(-dpi) / (tanh(-dpi) - tanh(dpi))
    stretch = cte_value / (tanh(dpi) - tanh(-dpi))

    if (sleft - 0.5_dp * width <= x .and. x <= sleft + 0.5_dp * width) then
      get_dropoff_dr = (2.0_dp * dpi * stretch / width) / cosh( &
        2.0_dp * dpi * (x - sleft) / width &
      )**2
    else if (sright - 0.5_dp * width <= x .and. x <= sright + 0.5_dp * width) then
      get_dropoff_dr = (-2.0_dp * dpi * stretch / width) / cosh( &
        2.0_dp * dpi * (sright - x) / width &
      )**2
    else
      get_dropoff_dr = 0.0_dp
    end if
  end function get_dropoff_dr

end module mod_physics_utils
