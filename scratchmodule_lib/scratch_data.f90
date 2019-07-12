module scratch_data
  use, intrinsic :: iso_c_binding
  implicit none
  private
  !alignment = 2kib
  !pad with 256 doubles == 2kib
  integer,parameter,public :: dp=selected_real_kind(14,200), &
                              int_long=selected_int_kind(16),&
                              alignment = 2048,&
                              sizeof_dp=c_sizeof(1.0_dp),&
                              alignment_in_real_dp = alignment/sizeof_dp,&
                              pad=256
  type :: scr
     real(dp), pointer,contiguous :: alignednode(:),node(:)
  end type scr
  type(scr),public,allocatable :: scratch(:),scratch_swap(:)
  type ::  reg
     integer(int_long) :: nodesize
     integer(int_long),pointer,contiguous :: usage(:,:)
     character(255),pointer,contiguous :: tag(:)
  end type reg
  type(reg),public,allocatable:: register(:),register_swap(:)
  !size(reg,1) .eq. size(node,1)
  !nodesize=size(scratch(id)%node)
  !size(usage,1)=> number of segments
  !size(usage,2)=5 => 1:size, 2:start, 3:end, 4:locked(0/1), 5:inuse(0/1)
  !size(tag,1)  => number of segments

end module scratch_data
