module scratch_segment_managment
  implicit none
  private

  public :: search_segment
  public :: init_finalize

contains
  subroutine init_finalize(task)
    use scratch_data, only: register
    use scratch_managment, only: init_scratch,finalize_scratch
    integer,intent(in) :: task
    integer :: i,n
    logical :: success
    if(task.eq.-1)then
       success=init_scratch()
       allocate(register(1))
       allocate(register(1)%usage(1,5))
       allocate(register(1)%tag(1))
       success=update_usage(1,1,locked=0,inuse=0,tag='',len=1,start=1,end=1)
       register(1)%nodesize=1
    elseif(task.eq.-2)then
       success=finalize_scratch()
       n=size(register)
       do i=1,n
          deallocate(register(i)%usage)
          deallocate(register(i)%tag)
       end do
       deallocate(register)
    end if
  end subroutine init_finalize

  subroutine search_segment(len,tag,task,arrayin,arrayout)
    use scratch_data, only : dp,register,pad,alignment_in_real_dp,int_long
    integer(int_long),intent(in) :: len
    integer,intent(in) :: task
    character(*),intent(in) :: tag
    real(dp),pointer,intent(in),optional,contiguous :: arrayin(:)
    real(dp),pointer,intent(out),optional,contiguous :: arrayout(:)
    integer :: id,id1,idsplit,id1split,num_nodes,num_segments,i,j
    integer(int_long) :: paddedlen
    logical :: split,success

    j=0
    i=0
    !pad memory
    paddedlen=len+pad
    !pad additionally so that we keep alignment
    paddedlen=paddedlen+mod(paddedlen,alignment_in_real_dp)
    num_nodes=size(register,1)
    if(task.eq.1)then
       split=.false.
       loop1_id: do id=1,num_nodes
          if (paddedlen.le.register(id)%nodesize)then
             !this node might hold big enough segments
             num_segments=size(register(id)%usage,1)
             loop1_id1: do id1=1,num_segments
                if (register(id)%usage(id1,4).eq.0)then
                   if(paddedlen.eq.register(id)%usage(id1,1))then
                      split=.false.
                      exit loop1_id
                   elseif(.not.split.and.paddedlen.lt.register(id)%usage(id1,1))then
                      split=.true.
                      idsplit=id
                      id1split=id1
                   end if
                end if
             end do loop1_id1
          end if
       end do loop1_id
       !found something usefull?
       if (split) then
          id=idsplit
          id1=id1split
          call split_segments(paddedlen,id,id1)
       !need new node?
       elseif(id.gt.num_nodes)then
          call combine_nodes(paddedlen,id,id1)
       end if
       !else: perfect match
       success=add_association(paddedlen,id,id1,arrayout)
       success=update_usage(id,id1,locked=1,inuse=1,tag=tag)

    elseif(task.eq.2.or.task.eq.3)then
       loop2_id: do id=1,num_nodes
          if (paddedlen.le.register(id)%nodesize)then
             !this node might be big enough
             num_segments=size(register(id)%usage,1)
             loop2_id1: do id1=1,num_segments
                if (register(id)%usage(id1,4).eq.1.and.&
                     register(id)%usage(id1,5).eq.1)then
                   if(paddedlen.eq.register(id)%usage(id1,1)&
                        .and.tag.eq.register(id)%tag(id1))then
                      !only exit if arrayin is associated with the found segment
                      !check association with request !!!len!!!
                      if (check_association(len,id,id1,arrayin)) exit loop2_id
                   end if
                end if
             end do loop2_id1
          end if
       end do loop2_id
       if(task.eq.2)then
          success=update_usage(id,id1,locked=0,inuse=0)
          call combine_segments(id,id1)
          call combine_nodes(0,id,id1)
       elseif(task.eq.3)then
          success=update_usage(id,id1,locked=1,inuse=0)
          call combine_nodes(0,id,id1)
       end if
    elseif(task.eq.4)then
       loop4_id: do id=1,num_nodes
          if (paddedlen.le.register(id)%nodesize)then
             !this node might be big enough
             num_segments=size(register(id)%usage,1)
             loop4_id1: do id1=1,num_segments
                if (register(id)%usage(id1,4).eq.1.and.&
                     register(id)%usage(id1,5).eq.0)then
                   !check also the tag of the possible seg
                   if(paddedlen.eq.register(id)%usage(id1,1)&
                        .and.tag.eq.register(id)%tag(id1))exit loop4_id
                end if
             end do loop4_id1
          end if
       end do loop4_id
       success=update_usage(id,id1,inuse=1,locked=1,tag=tag)
       success=add_association(paddedlen,id,id1,arrayout)
    end if
  end subroutine search_segment

  function check_association(len,id,id1,array) result(success)
    use scratch_data, only: register,dp,int_long
    use scratch_managment, only: add_remove_association
    logical ::success
    integer,intent(in):: id,id1
    integer(int_long), intent(in) :: len
    real(dp), pointer, intent(in),contiguous :: array(:)
    integer(int_long) :: start
    start=register(id)%usage(id1,2)
    call add_remove_association(len,id,start,2,success,arrayin=array)
  end function check_association

  function add_association(len,id,id1,array) result(success)
    use scratch_data, only: register,dp,int_long
    use scratch_managment, only: add_remove_association
    logical ::success
    integer(int_long),intent(in):: len
    integer,intent(in) :: id,id1
    real(dp), pointer, intent(out),contiguous :: array(:)
    integer(int_long) :: start
    start=register(id)%usage(id1,2)
    call add_remove_association(len,id,start,1,success,arrayout=array)
  end function add_association

  subroutine split_segments(newlen,id,id1)
    use scratch_data, only: register, int_long
    use scratch_utils, only: update_array
    !always split id1 into id1 and id1+1, new segment is always added as id1
    integer(int_long),intent(in) :: newlen
    integer,intent(in) :: id,id1
    integer(int_long),pointer,contiguous :: swapi(:,:)
    character(255),pointer,contiguous :: swapc(:)
    integer :: i,n
    integer(int_long) :: vali(5,2)
    character(255) :: valc(2)
    !add register entries
    n=size(register(id)%usage,1)+1
    allocate(swapi(n,5))
    allocate(swapc(n))
    !set new id1 settings
    vali(1,1)=newlen
    vali(2,1)=register(id)%usage(id1,2)
    vali(3,1)=vali(2,1)+vali(1,1)-1
    vali(4,1)=0
    vali(5,1)=0
    if(n.gt.1)then
       !set new id1+1 settings
       vali(1,2)=register(id)%usage(id1,1)-newlen
       vali(2,2)=vali(3,1)+1
       vali(3,2)=vali(2,2)+vali(1,2)-1
       vali(4,2)=0
       vali(5,2)=0
       valc=''
       call update_array(register(id)%usage,swapi,vali,id1)
       call update_array(register(id)%tag,swapc,valc,id1)
    else
       call update_array(register(id)%usage,swapi,vali(:,1),id1)
       call update_array(register(id)%tag,swapc,valc(1),id1)
    end if
    deallocate(register(id)%usage)
    deallocate(register(id)%tag)
    register(id)%usage=>swapi
    register(id)%tag=>swapc
  end subroutine split_segments

  subroutine combine_segments(id,id1)
    use scratch_data, only: register, int_long
    use scratch_utils, only: update_array
    !always combine id1 and id1+1 into id1
    integer,intent(in) :: id,id1
    integer(int_long),pointer,contiguous :: swapi(:,:)
    character(255),pointer,contiguous :: swapc(:)
    integer :: i,n
    integer(int_long) :: vali(5)
    character(255) :: valc
    !remove register entries
    n=size(register(id)%usage,1)
    !is id1 last segment?
    if(id1.eq.n.or.n.eq.1)return
    !id1+1 locked?
    if(register(id)%usage(id1+1,4).eq.1)return
    n=n-1
    allocate(swapi(n,5))
    allocate(swapc(n))
    !set new id1 settings
    vali(1)=register(id)%usage(id1,1)+register(id)%usage(id1+1,1)
    vali(2)=register(id)%usage(id1,2)
    vali(3)=vali(2)+vali(1)-1
    vali(4)=0
    vali(5)=0
    valc=''
    call update_array(register(id)%usage,swapi,vali,id1)
    call update_array(register(id)%tag,swapc,valc,id1)
    deallocate(register(id)%usage)
    deallocate(register(id)%tag)
    register(id)%usage=>swapi
    register(id)%tag=>swapc
  end subroutine combine_segments

  subroutine move_segment(src_id,src_id1,dest_id,dest_id1,len)
    use scratch_data, only: register,int_long
    use scratch_managment, only: move_scratch
    !move data from srcid srcid1 to destid destid1
    integer,intent(in) :: src_id,src_id1,dest_id,dest_id1
    integer(int_long), intent(in) :: len
    integer(int_long) :: src_start,dest_start
    logical :: success

    src_start=register(src_id)%usage(src_id1,2)
    dest_start=register(dest_id)%usage(dest_id1,2)
    success=move_scratch(src_id,dest_id,src_start,dest_start,len)
    success=update_usage(dest_id,dest_id1,locked=1,inuse=0,tag=register(src_id)%tag(src_id1))
    success=update_usage(src_id,src_id1,locked=0,inuse=0,tag='')
  end subroutine move_segment

  subroutine combine_nodes(len,id,id1)
    use scratch_data, only: register,int_long
    use scratch_managment, only: remove_scratch
    integer(int_long),intent(in) :: len
    integer,intent(out) :: id,id1
    integer(int_long) :: freed_len,newlen
    integer :: num_nodes
    logical :: success
    newlen=0
    num_nodes=size(register,1)
    !quick return?
    if(.not.(num_nodes.eq.1.and.register(1)%nodesize.ge.len))then
       !reverse loop
       do id=num_nodes,1,-1
          if(sum(register(id)%usage(:,4)).eq.0) then
             newlen=newlen+register(id)%nodesize
             success=remove_scratch(id)
             success=remove_register(id)
          end if
       end do
    end if
    if(len.gt.newlen)newlen=len
    call combine_free_moveable(newlen,len,id,id1)
  end subroutine combine_nodes

  subroutine combine_free_moveable(minlen,segmentlen,id_out,id1_out)
    use scratch_managment, only: add_scratch,remove_scratch
    use scratch_data, only: register, int_long
    integer,intent(out) :: id_out,id1_out
    integer(int_long),intent(in) :: minlen,segmentlen
    integer :: num_nodes,id,id1,count,num_segments
    integer(int_long) :: lenfree,lenmove,newlen,len
    logical,allocatable :: index_moveable(:)
    logical :: success
    lenfree=0
    lenmove=0
    num_nodes=size(register,1)
    allocate(index_moveable(num_nodes))
    index_moveable=.false.
    do id=1,num_nodes
       if(sum(register(id)%usage(:,5)).eq.0) then
          !this node is freely moveable,something must be locked since we cleared all
          !unlocked nodes bevore
          num_segments=size(register(id)%usage,1)
          do id1=1,num_segments
             if(register(id)%usage(id1,4).eq.1)then
                lenmove=lenmove+register(id)%usage(id1,1)
                index_moveable(id)=.true.
             elseif(register(id)%usage(id1,4).eq.0)then
                lenfree=lenfree+register(id)%usage(id1,1)
             end if
          end do
       end if
    end do
    !quick return?
    if(minlen.eq.0)then
       count=0
       do id=1,num_nodes
          if(index_moveable(id))count=count+1
       end do
       if(count.le.1)return
    end if
    !get a new scratch node
    if(minlen.gt.lenfree)then
       newlen=minlen+lenmove
    elseif(minlen.le.lenfree)then
       newlen=lenfree+lenmove
    end if
    success=add_scratch(newlen)
    success=add_register()
    register(1)%nodesize=newlen
    success=update_usage(1,1,locked=0,inuse=0,tag='',len=newlen,start=1,end=newlen)
    !do we need to move data?
    count=0
    if(lenmove.gt.0)then
       !reverse loop
       do id=num_nodes+1,2,-1
          !we added a new scratch node at position id=1, so old ids need to be incremented by 1
          if(index_moveable(id-1))then
             num_segments=size(register(id)%usage,1)
             do id1=1,num_segments
                if(register(id)%usage(id1,4).eq.1)then
                   len=register(id)%usage(id1,1)
                   count=count+1
                   call split_segments(len,1,count)
                   call move_segment(id,id1,1,count,len)
                end if
             end do
             success=remove_scratch(id)
             success=remove_register(id)
          end if
       end do
    end if
    count=count+1
    if(newlen-lenmove.gt.segmentlen.and.segmentlen.gt.0) call split_segments(segmentlen,1,count)
    id_out=1
    id1_out=count
  end subroutine combine_free_moveable

  function update_usage(id,id1,inuse,locked,tag,len,start,end) result(success)
    use scratch_data, only: register, int_long
    integer,intent(in) :: id,id1
    integer,intent(in),optional :: locked,inuse
    integer(int_long), intent(in), optional :: len,start,end
    character(*),intent(in),optional :: tag
    logical :: success
    success=.false.
    if(present(locked)) register(id)%usage(id1,4)=locked
    if(present(inuse)) register(id)%usage(id1,5)=inuse
    if(present(tag)) register(id)%tag(id1)=tag
    if(present(len)) register(id)%usage(id1,1)=len
    if(present(start)) register(id)%usage(id1,2)=start
    if(present(end)) register(id)%usage(id1,3)=end
    success=.true.
  end function update_usage

  function add_register() result(success)
    use scratch_data, only: register,register_swap
    logical :: success
    integer :: i,num_nodes
    success=.false.
    num_nodes=size(register,1)+1
    allocate(register_swap(num_nodes))
    allocate(register_swap(1)%usage(1,5))
    allocate(register_swap(1)%tag(1))
    do i=2,num_nodes
       register_swap(i)%nodesize = register(i-1)%nodesize
       register_swap(i)%usage => register(i-1)%usage
       register_swap(i)%tag => register(i-1)%tag
    end do
    deallocate(register)
    allocate(register(num_nodes))
    do i=1,num_nodes
       register(i)%nodesize = register_swap(i)%nodesize
       register(i)%usage => register_swap(i)%usage
       register(i)%tag => register_swap(i)%tag
    end do
    deallocate(register_swap)
    success=.true.
  end function add_register

  function remove_register(id) result(success)
    use scratch_data, only: register,register_swap
    integer,intent(in) :: id
    logical :: success
    integer :: i,num_nodes
    success=.false.
    num_nodes=size(register,1)-1
    allocate(register_swap(num_nodes))
    do i=1,id-1
       register_swap(i)%nodesize=register(i)%nodesize
       register_swap(i)%usage=> register(i)%usage
       register_swap(i)%tag=> register(i)%tag
    end do
    deallocate(register(id)%usage)
    deallocate(register(id)%tag)
    do i=id,num_nodes
       register_swap(i)%nodesize=register(i+1)%nodesize
       register_swap(i)%usage=> register(i+1)%usage
       register_swap(i)%tag=> register(i+1)%tag
    end do
    deallocate(register)
    allocate(register(num_nodes))
    do i=1,num_nodes
       register(i)%nodesize=register_swap(i)%nodesize
       register(i)%usage=>register_swap(i)%usage
       register(i)%tag=>register_swap(i)%tag
    end do
    deallocate(register_swap)
    success=.true.
  end function remove_register
!!!done
end module scratch_segment_managment

module scratch_segment_managment_interface
  use scratch_segment_managment, only: search_segment,init_finalize
  implicit none
  private
  public :: segment_interface

contains
  subroutine segment_interface(len,tag,task,arrayin,arrayout)
    use scratch_data, only : dp, int_long
    integer(int_long),intent(in),optional :: len
    integer,intent(in),optional :: task
    character(255),intent(in),optional :: tag
    real(dp),pointer,intent(in),optional,contiguous :: arrayin(:)
    real(dp),pointer,intent(out),optional,contiguous :: arrayout(:)
    logical,save :: first=.true.
    if(first) then
       call init_finalize(-1)
       first=.false.
    end if
    if(task.gt.0)then
       if(task.eq.1)then
          call search_segment(len,tag,task,arrayout=arrayout)
       elseif(task.eq.2.or.task.eq.3)then
          call search_segment(len,tag,task,arrayin=arrayin)
       elseif(task.eq.4)then
          call search_segment(len,tag,task,arrayout=arrayout)
       end if
    elseif(task.eq.-2)then
       call init_finalize(task)
    end if

  end subroutine segment_interface
end module scratch_segment_managment_interface
