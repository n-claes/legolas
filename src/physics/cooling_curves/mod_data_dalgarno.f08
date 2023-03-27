module mod_data_dalgarno
  use mod_global_variables, only: dp
  implicit none

  public

  integer, parameter :: n_dalgarno = 71
  !> log10 temperature values from Dalgarno and McCray (1978)
  real(dp), protected :: logT_dalgarno(n_dalgarno)
  !> log10 luminosity values from Dalgarno and McCray (1978)
  real(dp), protected :: logL_dalgarno(n_dalgarno)

  data logT_dalgarno / &
    2.0_dp, 2.1_dp, 2.2_dp, 2.3_dp, 2.4_dp, 2.5_dp, 2.6_dp, 2.7_dp, 2.8_dp, 2.9_dp, &
    3.0_dp, 3.1_dp, 3.2_dp, 3.3_dp, 3.4_dp, 3.5_dp, 3.6_dp, 3.7_dp, 3.8_dp, 3.9_dp, &
    4.0_dp, 4.1_dp, 4.2_dp, 4.3_dp, 4.4_dp, 4.5_dp, 4.6_dp, 4.7_dp, 4.8_dp, 4.9_dp, &
    5.0_dp, 5.1_dp, 5.2_dp, 5.3_dp, 5.4_dp, 5.5_dp, 5.6_dp, 5.7_dp, 5.8_dp, 5.9_dp, &
    6.0_dp, 6.1_dp, 6.2_dp, 6.3_dp, 6.4_dp, 6.5_dp, 6.6_dp, 6.7_dp, 6.8_dp, 6.9_dp, &
    7.0_dp, 7.1_dp, 7.2_dp, 7.3_dp, 7.4_dp, 7.5_dp, 7.6_dp, 7.7_dp, 7.8_dp, 7.9_dp, &
    8.0_dp, 8.1_dp, 8.2_dp, 8.3_dp, 8.4_dp, 8.5_dp, 8.6_dp, 8.7_dp, 8.8_dp, 8.9_dp, &
    9.0_dp &
  /
  data logL_dalgarno / &
    -26.523_dp, -26.398_dp, -26.301_dp, -26.222_dp, -26.097_dp, -26.011_dp, &
    -25.936_dp, -25.866_dp, -25.807_dp, -25.754_dp, -25.708_dp, -25.667_dp, &
    -25.630_dp, -25.595_dp, -25.564_dp, -25.534_dp, -25.506_dp, -25.479_dp, &
    -25.453_dp, -25.429_dp, -25.407_dp, -23.019_dp, -21.762_dp, -21.742_dp, &
    -21.754_dp, -21.730_dp, -21.523_dp, -21.455_dp, -21.314_dp, -21.229_dp, &
    -21.163_dp, -21.126_dp, -21.092_dp, -21.060_dp, -21.175_dp, -21.280_dp, &
    -21.390_dp, -21.547_dp, -21.762_dp, -22.050_dp, -22.271_dp, -22.521_dp, &
    -22.646_dp, -22.660_dp, -22.676_dp, -22.688_dp, -22.690_dp, -22.662_dp, &
    -22.635_dp, -22.609_dp, -22.616_dp, -22.646_dp, -22.697_dp, -22.740_dp, &
    -22.788_dp, -22.815_dp, -22.785_dp, -22.754_dp, -22.728_dp, -22.703_dp, &
    -22.680_dp, -22.630_dp, -22.580_dp, -22.530_dp, -22.480_dp, -22.430_dp, &
    -22.380_dp, -22.330_dp, -22.280_dp, -22.230_dp, -22.180_dp &
  /

end module mod_data_dalgarno
