module mod_cooling_curve_names
  use mod_global_variables, only: str_len
  implicit none

  public

  character(len=str_len), parameter :: JC_CORONA = "jc_corona"
  character(len=str_len), parameter :: DALGARNO = "dalgarno"
  character(len=str_len), parameter :: DALGARNO2 = "dalgarno2"
  character(len=str_len), parameter :: ML_SOLAR = "ml_solar"
  character(len=str_len), parameter :: SPEX = "spex"
  character(len=str_len), parameter :: SPEX_DALGARNO = "spex_dalgarno"
  character(len=str_len), parameter :: ROSNER = "rosner"
  character(len=str_len), parameter :: COLGAN = "colgan"
  character(len=str_len), parameter :: COLGAN_DM = "colgan_dm"

  character(len=str_len), parameter :: KNOWN_CURVES(9) = [ &
    JC_CORONA, DALGARNO, DALGARNO2, ML_SOLAR, SPEX, &
    SPEX_DALGARNO, ROSNER, COLGAN, COLGAN_DM &
  ]
end module mod_cooling_curve_names
