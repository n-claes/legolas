module mod_data_spex_enh
  use mod_global_variables, only: dp
  use mod_data_spex, only: n_spex
  implicit none

  public

  integer, parameter :: n_spex_enh_dalgarno = 76
  !> log10 temperature values to enhance the SPEX curve with Dalgarno and McCray (1978)
  real(dp), protected :: logT_spex_enh_dalgarno(n_spex_enh_dalgarno)
  !> log10 luminosity values to enhance the SPEX curve with Dalgarno and McCray (1978)
  real(dp), protected :: logL_spex_enh_dalgarno(n_spex_enh_dalgarno)
  !> luminosity values to enhance the SPEX curve
  real(dp), protected :: L_spex_enh(n_spex)

  data logT_spex_enh_dalgarno / &
    1.00_dp, 1.04_dp, 1.08_dp, 1.12_dp, 1.16_dp, 1.20_dp, 1.24_dp, 1.28_dp, 1.32_dp, &
    1.36_dp, 1.40_dp, 1.44_dp, 1.48_dp, 1.52_dp, 1.56_dp, 1.60_dp, 1.64_dp, 1.68_dp, &
    1.72_dp, 1.76_dp, 1.80_dp, 1.84_dp, 1.88_dp, 1.92_dp, 1.96_dp, 2.00_dp, 2.04_dp, &
    2.08_dp, 2.12_dp, 2.16_dp, 2.20_dp, 2.24_dp, 2.28_dp, 2.32_dp, 2.36_dp, 2.40_dp, &
    2.44_dp, 2.48_dp, 2.52_dp, 2.56_dp, 2.60_dp, 2.64_dp, 2.68_dp, 2.72_dp, 2.76_dp, &
    2.80_dp, 2.84_dp, 2.88_dp, 2.92_dp, 2.96_dp, 3.00_dp, 3.04_dp, 3.08_dp, 3.12_dp, &
    3.16_dp, 3.20_dp, 3.24_dp, 3.28_dp, 3.32_dp, 3.36_dp, 3.40_dp, 3.44_dp, 3.48_dp, &
    3.52_dp, 3.56_dp, 3.60_dp, 3.64_dp, 3.68_dp, 3.72_dp, 3.76_dp, 3.80_dp, 3.84_dp, &
    3.88_dp, 3.92_dp, 3.96_dp, 4.00_dp &
  /
  data logL_spex_enh_dalgarno / &
    -30.0377_dp, -29.7062_dp, -29.4055_dp, -29.1331_dp, -28.8864_dp, -28.6631_dp, &
    -28.4614_dp, -28.2791_dp, -28.1146_dp, -27.9662_dp, -27.8330_dp, -27.7129_dp, &
    -27.6052_dp, -27.5088_dp, -27.4225_dp, -27.3454_dp, -27.2767_dp, -27.2153_dp, &
    -27.1605_dp, -27.1111_dp, -27.0664_dp, -27.0251_dp, -26.9863_dp, -26.9488_dp, &
    -26.9119_dp, -26.8742_dp, -26.8353_dp, -26.7948_dp, -26.7523_dp, -26.7080_dp, &
    -26.6619_dp, -26.6146_dp, -26.5666_dp, -26.5183_dp, -26.4702_dp, -26.4229_dp, &
    -26.3765_dp, -26.3317_dp, -26.2886_dp, -26.2473_dp, -26.2078_dp, -26.1704_dp, &
    -26.1348_dp, -26.1012_dp, -26.0692_dp, -26.0389_dp, -26.0101_dp, -25.9825_dp, &
    -25.9566_dp, -25.9318_dp, -25.9083_dp, -25.8857_dp, -25.8645_dp, -25.8447_dp, &
    -25.8259_dp, -25.8085_dp, -25.7926_dp, -25.7778_dp, -25.7642_dp, -25.7520_dp, &
    -25.7409_dp, -25.7310_dp, -25.7222_dp, -25.7142_dp, -25.7071_dp, -25.7005_dp, &
    -25.6942_dp, -25.6878_dp, -25.6811_dp, -25.6733_dp, -25.6641_dp, -25.6525_dp, &
    -25.6325_dp, -25.6080_dp, -25.5367_dp, -25.4806_dp &
  /
  data L_spex_enh / &
    0.000013264_dp, 0.000042428_dp, 0.000088276_dp, 0.00017967_dp, &
    0.00084362_dp, 0.0034295_dp, 0.013283_dp, 0.042008_dp, &
    0.12138_dp, 0.30481_dp, 0.53386_dp, 0.76622_dp, &
    0.89459_dp, 0.95414_dp, 0.98342_dp, &
    1.0046_dp, 1.0291_dp, 1.0547_dp, 1.0767_dp, 1.0888_dp, 1.0945_dp, 1.0972_dp, &
    1.0988_dp, 1.1004_dp, 1.1034_dp, 1.1102_dp, 1.1233_dp, 1.1433_dp, 1.1638_dp, &
    1.1791_dp, 1.1885_dp, 1.1937_dp, 1.1966_dp, 1.1983_dp, 1.1993_dp, 1.1999_dp, &
    1.2004_dp, 1.2008_dp, 1.2012_dp, 1.2015_dp, 1.2020_dp, 1.2025_dp, 1.2030_dp, &
    1.2035_dp, 1.2037_dp, 1.2039_dp, 1.2040_dp, 1.2041_dp, 1.2042_dp, 1.2044_dp, &
    1.2045_dp, 1.2046_dp, 1.2047_dp, 1.2049_dp, 1.2050_dp, 1.2051_dp, 1.2053_dp, &
    1.2055_dp, 1.2056_dp, 1.2058_dp, 1.2060_dp, 1.2062_dp, 1.2065_dp, 1.2067_dp, &
    1.2070_dp, 1.2072_dp, 1.2075_dp, 1.2077_dp, 1.2078_dp, 1.2079_dp, 1.2080_dp, &
    1.2081_dp, 1.2082_dp, 1.2083_dp, 1.2083_dp, 1.2084_dp, 1.2084_dp, 1.2085_dp, &
    1.2085_dp, 1.2086_dp, 1.2086_dp, 1.2087_dp, 1.2087_dp, 1.2088_dp, 1.2088_dp, &
    1.2089_dp, 1.2089_dp, 1.2089_dp, 1.2089_dp, 1.2089_dp, 1.2090_dp, 1.2090_dp, &
    1.2090_dp, 1.2090_dp, 1.2090_dp, 1.2090_dp, 1.2090_dp, 1.2090_dp, 1.2090_dp, &
    1.2090_dp, 1.2090_dp, 1.2090_dp, 1.2090_dp, 1.2090_dp, 1.2090_dp, 1.2090_dp, &
    1.2090_dp, 1.2090_dp, 1.2090_dp, 1.2090_dp &
  /

end module mod_data_spex_enh
