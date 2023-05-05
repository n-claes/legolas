module mod_dims
  implicit none

  private

  type, public :: dims_t
    integer, private :: dim_integralblock
    integer, private :: dim_subblock
    integer, private :: dim_quadblock
    integer, private :: dim_matrix

  contains

    procedure, public :: set_block_dims
    procedure, public :: get_dim_integralblock
    procedure, public :: get_dim_subblock
    procedure, public :: get_dim_quadblock
    procedure, public :: get_dim_matrix

  end type dims_t

  public :: new_block_dims

contains

  pure function new_block_dims() result(dims)
    type(dims_t) :: dims

    dims%dim_integralblock = 2
    dims%dim_subblock = 0
    dims%dim_quadblock = 0
    dims%dim_matrix = 0
  end function new_block_dims


  pure subroutine set_block_dims(this, nb_eqs, gridpts)
    class(dims_t), intent(inout) :: this
    integer, intent(in) :: nb_eqs
    integer, intent(in) :: gridpts

    this%dim_subblock = nb_eqs * this%dim_integralblock
    this%dim_quadblock = this%dim_integralblock * this%dim_subblock
    this%dim_matrix = gridpts * this%dim_subblock
  end subroutine set_block_dims


  pure integer function get_dim_integralblock(this)
    class(dims_t), intent(in) :: this
    get_dim_integralblock = this%dim_integralblock
  end function get_dim_integralblock


  pure integer function get_dim_subblock(this)
    class(dims_t), intent(in) :: this
    get_dim_subblock = this%dim_subblock
  end function get_dim_subblock


  pure integer function get_dim_quadblock(this)
    class(dims_t), intent(in) :: this
    get_dim_quadblock = this%dim_quadblock
  end function get_dim_quadblock


  pure integer function get_dim_matrix(this)
    class(dims_t), intent(in) :: this
    get_dim_matrix = this%dim_matrix
  end function get_dim_matrix

end module mod_dims
