module scratch_utils
  use scratch_data, only: int_long
  implicit none
  private
  public :: update_array
  !alywas combine index and index+1 into index
  !always split index into index and index+1

  interface update_array
     module procedure split_array_i2
     module procedure combine_array_i2
     module procedure split_array_c1
     module procedure combine_array_c1
  end interface update_array

contains
  subroutine split_array_i2(arrayin,arrayout,val,index)
    integer,intent(in) :: index
    integer(int_long),intent(in),contiguous :: val(:,:),arrayin(:,:)
    integer(int_long),intent(out), contiguous :: arrayout(:,:)
    integer :: i,j,n1,n2
    n2=size(arrayout,2)
    n1=size(arrayout,1)
    do i=1,n2
       do j=1,index-1
          arrayout(j,i)=arrayin(j,i)
       end do
       arrayout(index,i)=val(i,1)
       arrayout(index+1,i)=val(i,2)
       do j=index+2,n1
          arrayout(j,i)=arrayin(j-1,i)
       end do
    end do
  end subroutine split_array_i2

  subroutine combine_array_i2(arrayin,arrayout,val,index)
    integer,intent(in) :: index
    integer(int_long),intent(in),contiguous :: val(:),arrayin(:,:)
    integer(int_long),intent(out),contiguous :: arrayout(:,:)
    integer :: i,j,n1,n2
    n2=size(arrayout,2)
    n1=size(arrayout,1)
    do i=1,n2
       do j=1,index-1
          arrayout(j,i)=arrayin(j,i)
       end do
       arrayout(index,i)=val(i)
       do j=index+1,n1
          arrayout(j,i)=arrayin(j+1,i)
       end do
    end do
  end subroutine combine_array_i2

  subroutine split_array_c1(arrayin,arrayout,val,index)
    integer,intent(in) :: index
    character(255),intent(in),contiguous :: arrayin(:),val(:)
    character(255),intent(out),contiguous :: arrayout(:)
    integer :: i,n
    n=size(arrayout,1)
    do i=1,index-1
       arrayout(i)=arrayin(i)
    end do
    arrayout(index)=val(1)
    arrayout(index+1)=val(2)
    do i=index+2,n
       arrayout(i)=arrayin(i-1)
    end do
  end subroutine split_array_c1

  subroutine combine_array_c1(arrayin,arrayout,val,index)
    integer,intent(in) :: index
    character(255),intent(in) :: val
    character(255),intent(in),contiguous :: arrayin(:)
    character(255),intent(out),contiguous :: arrayout(:)
    integer :: i,n
    n=size(arrayout,1)
    do i=1,index-1
       arrayout(i)=arrayin(i)
    end do
    arrayout(index)=val
    do i=index+1,n
       arrayout(i)=arrayin(i+1)
    end do
  end subroutine combine_array_c1
end module scratch_utils
