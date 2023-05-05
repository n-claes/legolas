module mod_data_mlsolar
  use mod_global_variables, only: dp
  implicit none

  public

  integer, parameter :: n_mlsolar = 71
  real(dp), protected :: logT_mlsolar(n_mlsolar)
  real(dp), protected :: logL_mlsolar(n_mlsolar)

  !> log10 temperature values from Melemma and Lundqvist (2002), solar metallicity
  data logT_mlsolar / &
    2.0_dp, 2.1_dp, 2.2_dp, 2.3_dp, 2.4_dp, 2.5_dp, 2.6_dp, 2.7_dp, 2.8_dp, 2.9_dp, &
    3.0_dp, 3.1_dp, 3.2_dp, 3.3_dp, 3.4_dp, 3.5_dp, 3.6_dp, 3.7_dp, 3.8_dp, 3.9_dp, &
    4.0_dp, 4.1_dp, 4.2_dp, 4.3_dp, 4.4_dp, 4.5_dp, 4.6_dp, 4.7_dp, 4.8_dp, 4.9_dp, &
    5.0_dp, 5.1_dp, 5.2_dp, 5.3_dp, 5.4_dp, 5.5_dp, 5.6_dp, 5.7_dp, 5.8_dp, 5.9_dp, &
    6.0_dp, 6.1_dp, 6.2_dp, 6.3_dp, 6.4_dp, 6.5_dp, 6.6_dp, 6.7_dp, 6.8_dp, 6.9_dp, &
    7.0_dp, 7.1_dp, 7.2_dp, 7.3_dp, 7.4_dp, 7.5_dp, 7.6_dp, 7.7_dp, 7.8_dp, 7.9_dp, &
    8.0_dp, 8.1_dp, 8.2_dp, 8.3_dp, 8.4_dp, 8.5_dp, 8.6_dp, 8.7_dp, 8.8_dp, 8.9_dp, &
    9.0_dp &
  /

  !> log10 luminosity values from Melemma and Lundqvist (2002), solar metallicity
  data logL_mlsolar / &
    -26.983_dp, -26.951_dp, -26.941_dp, -26.940_dp, -26.956_dp, -26.980_dp, &
    -27.011_dp, -27.052_dp, -27.097_dp, -27.145_dp, -27.195_dp, -27.235_dp, &
    -27.279_dp, -27.327_dp, -27.368_dp, -27.415_dp, -27.456_dp, -27.485_dp, &
    -27.468_dp, -27.223_dp, -25.823_dp, -23.501_dp, -22.162_dp, -22.084_dp, &
    -22.157_dp, -22.101_dp, -21.974_dp, -21.782_dp, -21.542_dp, -21.335_dp, &
    -21.251_dp, -21.275_dp, -21.236_dp, -21.173_dp, -21.167_dp, -21.407_dp, &
    -21.670_dp, -21.788_dp, -21.879_dp, -22.008_dp, -22.192_dp, -22.912_dp, &
    -22.918_dp, -22.887_dp, -22.929_dp, -23.023_dp, -23.094_dp, -23.117_dp, &
    -23.108_dp, -23.083_dp, -23.049_dp, -23.011_dp, -22.970_dp, -22.928_dp, &
    -22.885_dp, -22.842_dp, -22.798_dp, -22.754_dp, -22.709_dp, -22.665_dp, &
    -22.620_dp, -22.570_dp, -22.520_dp, -22.470_dp, -22.420_dp, -22.370_dp, &
    -22.320_dp, -22.270_dp, -22.220_dp, -22.170_dp, -22.120_dp &
  /

end module mod_data_mlsolar
