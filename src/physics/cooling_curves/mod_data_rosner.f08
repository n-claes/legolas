module mod_data_rosner
  use mod_global_variables, only: dp
  implicit none

  public

  !> @note
  !!     data for the piecewise analytical cooling curve by Rosner, Tucker and Vaiana
  !!     (1978), extended by Priest (1982).
  !! @endnote

  !> @note
  !!     Original values for xi in the rosner curve are given in SI,
  !!     these here are scaled to cgs.
  !!     \( 1 Wm^{-3} = 10^{-13} erg/(s*cm^3) \), and since these are in a log10 scale
  !!     that means a -13 difference.
  !! @endnote

  real(dp), protected :: logT_rosner(8)
  real(dp), protected :: logxi_rosner(9)
  real(dp), protected :: alpha_rosner(9)

  data logT_rosner / &
    3.89063_dp, 4.30195_dp, 4.575_dp, 4.9_dp, 5.4_dp, 5.77_dp, 6.315_dp, 7.60457_dp &
  /
  data logxi_rosner / &
    -69.900_dp, -48.307_dp, -21.850_dp, -31.000_dp, -21.200_dp, -10.400_dp, &
    -21.940_dp, -17.730_dp, -26.602_dp &
  /
  data alpha_rosner / &
    11.7_dp, 6.15_dp, 0.0_dp, 2.0_dp, 0.0_dp, -2.0_dp, 0.0_dp, -0.666666667_dp, 0.5_dp &
  /

end module mod_data_rosner
