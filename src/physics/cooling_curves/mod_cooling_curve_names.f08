module mod_cooling_curve_names
  use mod_global_variables, only: str_len
  implicit none

  public

  character(len=str_len), parameter :: JC_CORONA = "jc_corona"
  character(len=str_len), parameter :: DALGARNO = "dalgarno"
  character(len=str_len), parameter :: ML_SOLAR = "ml_solar"
  character(len=str_len), parameter :: SPEX = "spex"
  character(len=str_len), parameter :: SPEX_DALGARNO = "spex_dalgarno"
  character(len=str_len), parameter :: ROSNER = "rosner"

  character(len=str_len), parameter :: KNOWN_CURVES(6) = [ &
    JC_CORONA, DALGARNO, ML_SOLAR, SPEX, SPEX_DALGARNO, ROSNER &
  ]
end module mod_cooling_curve_names
