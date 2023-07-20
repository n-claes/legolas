module mod_state_vector_names
  use mod_global_variables, only: str_len_arr
  implicit none

  public

  character(len=str_len_arr), parameter, public :: sv_RHO_name = "rho"
  character(len=str_len_arr), parameter, public :: sv_V1_name = "v1"
  character(len=str_len_arr), parameter, public :: sv_V2_name = "v2"
  character(len=str_len_arr), parameter, public :: sv_V3_name = "v3"
  character(len=str_len_arr), parameter, public :: sv_T_name = "T"
  character(len=str_len_arr), parameter, public :: sv_A1_name = "a1"
  character(len=str_len_arr), parameter, public :: sv_A2_name = "a2"
  character(len=str_len_arr), parameter, public :: sv_A3_name = "a3"

end module mod_state_vector_names
