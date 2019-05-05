subroutine mp_min_in_place_int4_r0(data,comm)
#ifdef __PARALLEL
  use mpi
#endif
  implicit none
  ! Arguments
  integer, intent(in) :: comm
  integer(int_4), intent(inout) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_min_in_place_int4_r0'
#ifdef __PARALLEL
  integer :: ierr
  integer(int_4) :: stack_buff
  call mpi_allreduce(data,stack_buff,1,MPI_INTEGER,&
         MPI_MIN,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  data=stack_buff
#endif
end subroutine mp_min_in_place_int4_r0
subroutine mp_min_in_place_int8_r0(data,comm)
#ifdef __PARALLEL
  use mpi
#endif
  implicit none
  ! Arguments
  integer, intent(in) :: comm
  integer(int_8), intent(inout) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_min_in_place_int8_r0'
#ifdef __PARALLEL
  integer :: ierr
  integer(int_8) :: stack_buff
  call mpi_allreduce(data,stack_buff,1,MPI_INTEGER8,&
         MPI_MIN,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  data=stack_buff
#endif
end subroutine mp_min_in_place_int8_r0
subroutine mp_min_in_place_real8_r0(data,comm)
#ifdef __PARALLEL
  use mpi
#endif
  implicit none
  ! Arguments
  integer, intent(in) :: comm
  real(real_8), intent(inout) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_min_in_place_real8_r0'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: stack_buff
  call mpi_allreduce(data,stack_buff,1,MPI_DOUBLE_PRECISION,&
         MPI_MIN,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  data=stack_buff
#endif
end subroutine mp_min_in_place_real8_r0
subroutine mp_min_in_place_complex8_r0(data,comm)
#ifdef __PARALLEL
  use mpi
#endif
  implicit none
  ! Arguments
  integer, intent(in) :: comm
  complex(real_8), intent(inout) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_min_in_place_complex8_r0'
#ifdef __PARALLEL
  integer :: ierr
  complex(real_8) :: stack_buff
  call mpi_allreduce(data,stack_buff,1,MPI_DOUBLE_COMPLEX,&
         MPI_MIN,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  data=stack_buff
#endif
end subroutine mp_min_in_place_complex8_r0
