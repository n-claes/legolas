! =============================================================================
!> All data of the different cooling curves is contained in this module,
!! along with handling piecewise cooling curves.
!! @note The currently implemented cooling curves are
!!
!! - _"jc_corona"_: Colgan et al. (2008)
!! - _"dalgarno"_: Dalgarno and McCray (1978)
!! - _"ml_solar"_: Mellema and Lundqvist (2002), solar metallicity
!! - _"spex"_: Schure et al. (2009)
!! - _"spex_dalgarno"_: spex curve extended by the dalgarno one at low temperatures.
!! - _"rosner"_: piecewise analytical cooling curve by Rosner, Tucker and Vaiana (1978),
!!               extended by Priest (1982)  @endnote
module mod_cooling_curves
  use mod_global_variables, only: dp
  implicit none

  !> amount of points in the jc_corona curve
  integer, protected  :: n_jc_corona
  !> amount of points in the dalgarno curve
  integer, protected  :: n_dalgarno
  !> amount of points in the ml_solar curve
  integer, protected  :: n_ml_solar
  !> amount of points in the spex curve
  integer, protected  :: n_spex
  !> amount of points to enhance the spex_dalgarno curve
  integer, protected  :: n_dalgarno2
  !> log temperature datapoints in the jc_corona curve
  real(dp), protected :: t_jc_corona(45)
  !> log temperature datapoints in the dalgarno curve
  real(dp), protected :: t_dalgarno(71)
  !> log temperature datapoints in the ml_solar curve
  real(dp), protected :: t_ml_solar(71)
  !> log temperature datapoints in the spex curve
  real(dp), protected :: t_spex(110)
  !> log temperature datapoints to enhance the spex_dalgarno curve
  real(dp), protected :: t_dalgarno2(76)
  !> log luminosity datapoints in the jc_corona curve
  real(dp), protected :: l_jc_corona(45)
  !> log luminosity datapoints in the dalgarno curve
  real(dp), protected :: l_dalgarno(71)
  !> log luminosity datapoints in the ml_solar curve
  real(dp), protected :: l_ml_solar(71)
  !> log luminosity datapoints in the spex curve
  real(dp), protected :: l_spex(110)
  !> log luminosity datapoints to enhance the spex_dalgarno curve
  real(dp), protected :: l_dalgarno2(76)
  !> log luminosity datapoints in the spex_dalgarno curve
  real(dp), protected :: n_spex_enh(110)

  !> datapoints for the piecewise rosner curve
  real(dp), protected :: rosner_log10_T(8)
  !> log xi values for the rosner cooling curve
  real(dp), protected :: rosner_log10_xi(9)
  !> alpha values for the rosner cooling curve
  real(dp), protected :: rosner_alpha(9)

  data    n_jc_corona / 45 /
  data    t_jc_corona / 4.00000, 4.14230, 4.21995, 4.29761, 4.37528, &
                        4.45294, 4.53061, 4.60827, 4.68593, 4.76359, &
                        4.79705, 4.83049, 4.86394, 4.89739, 4.93084, &
                        4.96428, 4.99773, 5.03117, 5.06461, 5.17574, &
                        5.28684, 5.39796, 5.50907, 5.62018, 5.73129, &
                        5.84240, 5.95351, 6.06461, 6.17574, 6.28684, &
                        6.39796, 6.50907, 6.62018, 6.73129, 6.84240, &
                        6.95351, 7.06461, 7.17574, 7.28684, 7.39796, &
                        7.50907, 7.62018, 7.73129, 7.84240, 7.95351  /
  data    l_jc_corona / -200.18883, -100.78630, -30.60384, -22.68481, -21.76445, &
                        -21.67936, -21.54218, -21.37958, -21.25172, -21.17584, &
                        -21.15783, -21.14491, -21.13527, -21.12837, -21.12485, &
                        -21.12439, -21.12642, -21.12802, -21.12548, -21.08965, &
                        -21.08812, -21.19542, -21.34582, -21.34839, -21.31701, &
                        -21.29072, -21.28900, -21.34104, -21.43122, -21.62448, &
                        -21.86694, -22.02897, -22.08051, -22.06057, -22.01973, &
                        -22.00000, -22.05161, -22.22175, -22.41452, -22.52581, &
                        -22.56914, -22.57486, -22.56151, -22.53969, -22.51490  /

  data    n_dalgarno / 71 /
  data    t_dalgarno / 2.0, 2.1, 2.2, 2.3, 2.4, &
                       2.5, 2.6, 2.7, 2.8, 2.9, &
                       3.0, 3.1, 3.2, 3.3, 3.4, &
                       3.5, 3.6, 3.7, 3.8, 3.9, &
                       4.0, 4.1, 4.2, 4.3, 4.4, &
                       4.5, 4.6, 4.7, 4.8, 4.9, &
                       5.0, 5.1, 5.2, 5.3, 5.4, &
                       5.5, 5.6, 5.7, 5.8, 5.9, &
                       6.0, 6.1, 6.2, 6.3, 6.4, &
                       6.5, 6.6, 6.7, 6.8, 6.9, &
                       7.0, 7.1, 7.2, 7.3, 7.4, &
                       7.5, 7.6, 7.7, 7.8, 7.9, &
                       8.0, 8.1, 8.2, 8.3, 8.4, &
                       8.5, 8.6, 8.7, 8.8, 8.9, &
                       9.0 /
  data    l_dalgarno / -26.523, -26.398, -26.301, -26.222, -26.097, &
                       -26.011, -25.936, -25.866, -25.807, -25.754, &
                       -25.708, -25.667, -25.630, -25.595, -25.564, &
                       -25.534, -25.506, -25.479, -25.453, -25.429, &
                       -25.407, -23.019, -21.762, -21.742, -21.754, &
                       -21.730, -21.523, -21.455, -21.314, -21.229, &
                       -21.163, -21.126, -21.092, -21.060, -21.175, &
                       -21.280, -21.390, -21.547, -21.762, -22.050, &
                       -22.271, -22.521, -22.646, -22.660, -22.676, &
                       -22.688, -22.690, -22.662, -22.635, -22.609, &
                       -22.616, -22.646, -22.697, -22.740, -22.788, &
                       -22.815, -22.785, -22.754, -22.728, -22.703, &
                       -22.680, -22.630, -22.580, -22.530, -22.480, &
                       -22.430, -22.380, -22.330, -22.280, -22.230, &
                       -22.180 /

  data    n_ml_solar / 71 /
  data    t_ml_solar / 2.0, 2.1, 2.2, 2.3, 2.4, &
                       2.5, 2.6, 2.7, 2.8, 2.9, &
                       3.0, 3.1, 3.2, 3.3, 3.4, &
                       3.5, 3.6, 3.7, 3.8, 3.9, &
                       4.0, 4.1, 4.2, 4.3, 4.4, &
                       4.5, 4.6, 4.7, 4.8, 4.9, &
                       5.0, 5.1, 5.2, 5.3, 5.4, &
                       5.5, 5.6, 5.7, 5.8, 5.9, &
                       6.0, 6.1, 6.2, 6.3, 6.4, &
                       6.5, 6.6, 6.7, 6.8, 6.9, &
                       7.0, 7.1, 7.2, 7.3, 7.4, &
                       7.5, 7.6, 7.7, 7.8, 7.9, &
                       8.0, 8.1, 8.2, 8.3, 8.4, &
                       8.5, 8.6, 8.7, 8.8, 8.9, &
                       9.0 /
  data    l_ml_solar / -26.983, -26.951, -26.941, -26.940, -26.956, &
                       -26.980, -27.011, -27.052, -27.097, -27.145, &
                       -27.195, -27.235, -27.279, -27.327, -27.368, &
                       -27.415, -27.456, -27.485, -27.468, -27.223, &
                       -25.823, -23.501, -22.162, -22.084, -22.157, &
                       -22.101, -21.974, -21.782, -21.542, -21.335, &
                       -21.251, -21.275, -21.236, -21.173, -21.167, &
                       -21.407, -21.670, -21.788, -21.879, -22.008, &
                       -22.192, -22.912, -22.918, -22.887, -22.929, &
                       -23.023, -23.094, -23.117, -23.108, -23.083, &
                       -23.049, -23.011, -22.970, -22.928, -22.885, &
                       -22.842, -22.798, -22.754, -22.709, -22.665, &
                       -22.620, -22.570, -22.520, -22.470, -22.420, &
                       -22.370, -22.320, -22.270, -22.220, -22.170, &
                       -22.120 /

  data    n_spex / 110 /
  data    t_spex / 3.80, 3.84, 3.88, 3.92, 3.96, &
                   4.00, 4.04, 4.08, 4.12, 4.16, &
                   4.20, 4.24, 4.28, 4.32, 4.36, &
                   4.40, 4.44, 4.48, 4.52, 4.56, &
                   4.60, 4.64, 4.68, 4.72, 4.76, &
                   4.80, 4.84, 4.88, 4.92, 4.96, &
                   5.00, 5.04, 5.08, 5.12, 5.16, &
                   5.20, 5.24, 5.28, 5.32, 5.36, &
                   5.40, 5.44, 5.48, 5.52, 5.56, &
                   5.60, 5.64, 5.68, 5.72, 5.76, &
                   5.80, 5.84, 5.88, 5.92, 5.96, &
                   6.00, 6.04, 6.08, 6.12, 6.16, &
                   6.20, 6.24, 6.28, 6.32, 6.36, &
                   6.40, 6.44, 6.48, 6.52, 6.56, &
                   6.60, 6.64, 6.68, 6.72, 6.76, &
                   6.80, 6.84, 6.88, 6.92, 6.96, &
                   7.00, 7.04, 7.08, 7.12, 7.16, &
                   7.20, 7.24, 7.28, 7.32, 7.36, &
                   7.40, 7.44, 7.48, 7.52, 7.56, &
                   7.60, 7.64, 7.68, 7.72, 7.76, &
                   7.80, 7.84, 7.88, 7.92, 7.96, &
                   8.00, 8.04, 8.08, 8.12, 8.16  /
  data    l_spex / -25.7331, -25.0383, -24.4059, -23.8288, -23.3027, &
                   -22.8242, -22.3917, -22.0067, -21.6818, -21.4529, &
                   -21.3246, -21.3459, -21.4305, -21.5293, -21.6138, &
                   -21.6615, -21.6551, -21.5919, -21.5092, -21.4124, &
                   -21.3085, -21.2047, -21.1067, -21.0194, -20.9413, &
                   -20.8735, -20.8205, -20.7805, -20.7547, -20.7455, &
                   -20.7565, -20.7820, -20.8008, -20.7994, -20.7847, &
                   -20.7687, -20.7590, -20.7544, -20.7505, -20.7545, &
                   -20.7888, -20.8832, -21.0450, -21.2286, -21.3737, &
                   -21.4573, -21.4935, -21.5098, -21.5345, -21.5863, &
                   -21.6548, -21.7108, -21.7424, -21.7576, -21.7696, &
                   -21.7883, -21.8115, -21.8303, -21.8419, -21.8514, &
                   -21.8690, -21.9057, -21.9690, -22.0554, -22.1488, &
                   -22.2355, -22.3084, -22.3641, -22.4033, -22.4282, &
                   -22.4408, -22.4443, -22.4411, -22.4334, -22.4242, &
                   -22.4164, -22.4134, -22.4168, -22.4267, -22.4418, &
                   -22.4603, -22.4830, -22.5112, -22.5449, -22.5819, &
                   -22.6177, -22.6483, -22.6719, -22.6883, -22.6985, &
                   -22.7032, -22.7037, -22.7008, -22.6950, -22.6869, &
                   -22.6769, -22.6655, -22.6531, -22.6397, -22.6258, &
                   -22.6111, -22.5964, -22.5816, -22.5668, -22.5519, &
                   -22.5367, -22.5216, -22.5062, -22.4912, -22.4753 /
  data    n_spex_enh / 0.000013264, 0.000042428, 0.000088276, 0.00017967, &
                       0.00084362, 0.0034295, 0.013283, 0.042008, &
                       0.12138, 0.30481, 0.53386, 0.76622, &
                       0.89459, 0.95414, 0.98342, &
                       1.0046, 1.0291, 1.0547, 1.0767, 1.0888, 1.0945, &
                       1.0972, 1.0988, 1.1004, 1.1034, 1.1102, 1.1233, &
                       1.1433, 1.1638, 1.1791, 1.1885, 1.1937, 1.1966, &
                       1.1983, 1.1993, 1.1999, 1.2004, 1.2008, 1.2012, &
                       1.2015, 1.2020, 1.2025, 1.2030, 1.2035, 1.2037, &
                       1.2039, 1.2040, 1.2041, 1.2042, 1.2044, 1.2045, &
                       1.2046, 1.2047, 1.2049, 1.2050, 1.2051, 1.2053, &
                       1.2055, 1.2056, 1.2058, 1.2060, 1.2062, 1.2065, &
                       1.2067, 1.2070, 1.2072, 1.2075, 1.2077, 1.2078, &
                       1.2079, 1.2080, 1.2081, 1.2082, 1.2083, 1.2083, &
                       1.2084, 1.2084, 1.2085, 1.2085, 1.2086, 1.2086, &
                       1.2087, 1.2087, 1.2088, 1.2088, 1.2089, 1.2089, &
                       1.2089, 1.2089, 1.2089, 1.2090, 1.2090, 1.2090, &
                       1.2090, 1.2090, 1.2090, 1.2090, 1.2090, 1.2090, &
                       1.2090, 1.2090, 1.2090, 1.2090, 1.2090, 1.2090, &
                       1.2090, 1.2090, 1.2090, 1.2090, 1.2090    /

  data    n_dalgarno2 / 76 /
  data    t_dalgarno2 / 1.00, 1.04, 1.08, 1.12, 1.16, &
                        1.20, 1.24, 1.28, 1.32, 1.36, &
                        1.40, 1.44, 1.48, 1.52, 1.56, &
                        1.60, 1.64, 1.68, 1.72, 1.76, &
                        1.80, 1.84, 1.88, 1.92, 1.96, &
                        2.00, 2.04, 2.08, 2.12, 2.16, &
                        2.20, 2.24, 2.28, 2.32, 2.36, &
                        2.40, 2.44, 2.48, 2.52, 2.56, &
                        2.60, 2.64, 2.68, 2.72, 2.76, &
                        2.80, 2.84, 2.88, 2.92, 2.96, &
                        3.00, 3.04, 3.08, 3.12, 3.16, &
                        3.20, 3.24, 3.28, 3.32, 3.36, &
                        3.40, 3.44, 3.48, 3.52, 3.56, &
                        3.60, 3.64, 3.68, 3.72, 3.76, &
                        3.80, 3.84, 3.88, 3.92, 3.96, &
                        4.00 /
  !! @note The spex_dalgarno curve assumes an ionisation fraction of <tt>1e-3</tt>. @endnote
  data    l_dalgarno2 / -30.0377, -29.7062, -29.4055, -29.1331, &
                        -28.8864, -28.6631, -28.4614, -28.2791, &
                        -28.1146, -27.9662, -27.8330, -27.7129, &
                        -27.6052, -27.5088, -27.4225, -27.3454, &
                        -27.2767, -27.2153, -27.1605, -27.1111, &
                        -27.0664, -27.0251, -26.9863, -26.9488, &
                        -26.9119, -26.8742, -26.8353, -26.7948, &
                        -26.7523, -26.7080, -26.6619, -26.6146, &
                        -26.5666, -26.5183, -26.4702, -26.4229, &
                        -26.3765, -26.3317, -26.2886, -26.2473, &
                        -26.2078, -26.1704, -26.1348, -26.1012, &
                        -26.0692, -26.0389, -26.0101, -25.9825, &
                        -25.9566, -25.9318, -25.9083, -25.8857, &
                        -25.8645, -25.8447, -25.8259, -25.8085, &
                        -25.7926, -25.7778, -25.7642, -25.7520, &
                        -25.7409, -25.7310, -25.7222, -25.7142, &
                        -25.7071, -25.7005, -25.6942, -25.6878, &
                        -25.6811, -25.6733, -25.6641, -25.6525, &
                        -25.6325, -25.6080, -25.5367, -25.4806  /

  data    rosner_log10_T / 3.89063d0, 4.30195d0, 4.575d0, 4.9d0, 5.4d0, 5.77d0, &
                           6.315d0, 7.60457d0 /
  !! @note Original values for xi in the rosner curve are given in SI, these here are scaled to cgs.
  !!       \( 1 Wm^{-3} = 10^{-13} erg/(s*cm^3) \), and since these are in a log10 scale
  !!       that means a -13 difference. @endnote
  data    rosner_log10_xi / -69.900d0, -48.307d0, -21.850d0, -31.000d0, -21.200d0, &
                           -10.400d0, -21.940d0, -17.730d0, -26.602d0           /
  !
  data    rosner_alpha / 11.7d0, 6.15d0, 0.0d0, 2.0d0, 0.0d0, -2.0d0, 0.0d0, -0.666666667d0, 0.5d0 /

contains


  !> Uses the piecewise rosner cooling curve to calculate the radiative cooling values
  !! based on the equilibrium temperature.
  !! @note    The temperature array T0 has to be normalised when calling this routine. @endnote
  !! @note    Both <tt>lambda_T</tt> and <tt>d_lambda_dT</tt> are normalised on exit. @endnote
  !! @warning To ensure a correct calculation the temperature normalisation must be set
  !!          to the desired value. @endwarning
  subroutine get_rosner_cooling(settings, T0, lambda_T, dlambda_dT)
    use mod_settings, only: settings_t

    type(settings_t), intent(in) :: settings
    !> the equilibrium temperature, normalised
    real(dp), intent(in)  :: T0(:)
    !> the cooling values corresponding to T0
    real(dp), intent(out) :: lambda_T(:)
    !> the derivative of the cooling values corresponding to T0
    real(dp), intent(out) :: dlambda_dT(:)
    real(dp) :: unit_temperature, unit_lambdaT, unit_dlambdaT_dT
    real(dp)    :: log10_T0
    integer     :: i, j, idx

    unit_temperature = settings%units%get_unit_temperature()
    unit_lambdaT = settings%units%get_unit_lambdaT()
    unit_dlambdaT_dT = unit_lambdaT / unit_temperature

    do i = 1, size(T0)
      idx = 1
      log10_T0 = dlog10(T0(i) * unit_temperature)

      ! look up value in tables
      if (log10_T0 > rosner_log10_T(8)) then
        idx = 9
      else
        do j = 1, size(rosner_log10_T)
          if (log10_T0 < rosner_log10_T(j)) then
            idx = j
            exit
          end if
        end do
      end if

      ! lambda_T = xi * T**alpha, so log(lambda_T) = log(xi) + alpha*log(T)
      lambda_T(i) = 10.0d0 ** (rosner_log10_xi(idx) + rosner_alpha(idx) * log10_T0)
      ! dlambda_dT = alpha * xi * T**(alpha - 1) hence
      !            = alpha * 10**(log(xi) + (alpha-1)*log(T))
      dlambda_dT(i) = rosner_alpha(idx) * 10.0d0 ** ( &
        rosner_log10_xi(idx) + (rosner_alpha(idx)-1) * log10_T0 &
      )
    end do

    ! renormalise variables
    lambda_T = lambda_T / unit_lambdaT
    dlambda_dT = dlambda_dT / unit_dlambdaT_dT
  end subroutine get_rosner_cooling

end module mod_cooling_curves
