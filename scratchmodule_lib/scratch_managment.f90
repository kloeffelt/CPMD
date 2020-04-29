module scratch_managment
  implicit none
  private

  public :: add_scratch
  public :: remove_scratch
  public :: add_remove_association
  public :: init_scratch
  public :: finalize_scratch
  public :: move_scratch

contains

  function init_scratch() result(success)
    use scratch_data, only : scratch
    logical :: success
    success=.false.
    allocate(scratch(1))
    allocate(scratch(1)%node(0))
    success=.true.
  end function init_scratch

  function finalize_scratch() result(success)
    use scratch_data, only : scratch
    integer :: i,n
    logical :: success
    success=.false.
    n=size(scratch)
    do i=1,n
       deallocate(scratch(i)%node)
    end do
    deallocate(scratch)
    success=.true.
  end function finalize_scratch

  function calc_shift(arrayin) result(shift)
    use,intrinsic :: iso_c_binding
    use scratch_data, only : dp,int_long,alignment,sizeof_dp,alignment_in_real_dp
    real(dp),intent(in),pointer,contiguous :: arrayin(:)
    type(c_ptr) :: c_addr
    integer(int_long) :: faddr, shift
    c_addr=c_loc(arrayin(1))
    faddr=transfer(c_addr,faddr)
    shift=(alignment-mod(faddr,alignment))/sizeof_dp
    if(shift.eq.alignment_in_real_dp) shift=0
  end function calc_shift

  function add_scratch(len) result(success)
    use scratch_data, only : scratch, scratch_swap, alignment_in_real_dp, int_long
    integer(int_long),intent(in) :: len
    integer(int_long):: shift
    integer :: i,n
    logical :: success
    success=.false.
    n=size(scratch)+1
    allocate(scratch_swap(n))
    allocate(scratch_swap(1)%node(len+alignment_in_real_dp))
    shift=calc_shift(scratch_swap(1)%node(:))
    scratch_swap(1)%alignednode(1:len)=> scratch_swap(1)%node(1+shift:len+shift)
    do i=2,n
       scratch_swap(i)%node=>scratch(i-1)%node
       scratch_swap(i)%alignednode=>scratch(i-1)%alignednode
    end do
    deallocate(scratch)
    allocate(scratch(n))
    do i=1,n
       scratch(i)%alignednode=>scratch_swap(i)%alignednode
       scratch(i)%node=>scratch_swap(i)%node
    end do
    deallocate(scratch_swap)
    success=.true.
  end function add_scratch

  function remove_scratch(id) result(success)
    use scratch_data, only : scratch, scratch_swap, int_long
    integer,intent(in) :: id
    integer(int_long) :: n, i
    logical :: success
    success=.false.
    n=size(scratch)-1
    allocate(scratch_swap(n))
    do i=1,id-1
       scratch_swap(i)%node=>scratch(i)%node
       scratch_swap(i)%alignednode=>scratch(i)%alignednode
    end do
    do i=id+1,n+1
       scratch_swap(i-1)%node=>scratch(i)%node
       scratch_swap(i-1)%alignednode=>scratch(i)%alignednode
    end do
    deallocate(scratch(id)%node)
    deallocate(scratch)
    allocate(scratch(n))
    do i=1,n
       scratch(i)%node=>scratch_swap(i)%node
       scratch(i)%alignednode=>scratch_swap(i)%alignednode
    end do
    deallocate(scratch_swap)
    success=.true.
  end function  remove_scratch

  subroutine add_remove_association(len,id,start,task,success,arrayin,arrayout)
    use scratch_data, only : scratch,dp,int_long
    integer,intent(in) :: id,task
    integer(int_long),intent(in) :: len, start
    logical :: success
    real(dp), pointer, intent(in), optional, contiguous :: arrayin(:)
    real(dp), pointer, intent(out), optional, contiguous :: arrayout(:)
    real(dp), pointer, contiguous :: temp(:)
    success=.false.
    if(task.eq.1.or.task.eq.4)then
       if(.not.present(arrayout)) stop "add_remove_association: missing arrayout"
       arrayout(1:len)=>scratch(id)%alignednode(start:start+len-1)
       success=.true.
    elseif(task.eq.2.or.task.eq.3)then
       if(.not.present(arrayin)) stop "add_remove_association: missing arrayin:"
       temp(1:len)=>scratch(id)%alignednode(start:start+len-1)
       if(associated(arrayin,temp)) success=.true.
    end if
  end subroutine add_remove_association

  function move_scratch(src_id,dest_id,src_start,dest_start,len) result(success)
    use scratch_data, only : scratch, int_long
    integer,intent(in) :: src_id,dest_id
    integer(int_long),intent(in) :: src_start,dest_start,len
    logical :: success
    integer(int_long) :: deststart,srcstart
    success=.false.
    deststart=dest_start-1
    srcstart=src_start-1
    !get rid of pointer attribute
    call copy_array(scratch(src_id)%alignednode,scratch(dest_id)%alignednode,srcstart,deststart,len,success)
  end function move_scratch

  pure subroutine copy_array(arrayin,arrayout,instart,outstart,len,success)
    use scratch_data, only : dp, int_long
    integer(int_long),intent(in) :: instart,outstart,len
    real(dp),contiguous,intent(in) :: arrayin(:)
    real(dp),contiguous,intent(out) :: arrayout(:)
    logical,intent(out) :: success
    integer(int_long) :: i
    success=.false.
    do i=1,len
       arrayout(outstart+i)=arrayin(instart+i)
    end do
    success=.true.
  end subroutine  copy_array

end module scratch_managment
