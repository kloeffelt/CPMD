subroutine mp_sum_root_int4_r0(data_in,data_out, root, comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) ::  root, comm
  integer(int_4), intent(in)  :: data_in
  integer(int_4), intent(inout)  :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_int4_r0'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+1*mp_int4_in_bytes
  tim1=m_walltime()      
  call mpi_reduce(data_in,data_out,1,MPI_INTEGER,&
       MPI_SUM,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out =data_in 
#endif
end subroutine mp_sum_root_int4_r0
subroutine mp_sum_root_int4_r1(data_in,data_out,n, root, comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root, comm
  integer(int_4), intent(in) , dimension(*) :: data_in
  integer(int_4), intent(inout) , dimension(*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_int4_r1'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_int4_in_bytes
  tim1=m_walltime()      
  call mpi_reduce(data_in,data_out,n,MPI_INTEGER,&
       MPI_SUM,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out (1:n)=data_in (1:n)
#endif
end subroutine mp_sum_root_int4_r1
subroutine mp_sum_root_int4_r2(data_in,data_out,n, root, comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root, comm
  integer(int_4), intent(in) , dimension(1,*) :: data_in
  integer(int_4), intent(inout) , dimension(1,*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_int4_r2'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_int4_in_bytes
  tim1=m_walltime()      
  call mpi_reduce(data_in,data_out,n,MPI_INTEGER,&
       MPI_SUM,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out (1,1:n)=data_in (1,1:n)
#endif
end subroutine mp_sum_root_int4_r2
subroutine mp_sum_root_real8_r0(data_in,data_out, root, comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) ::  root, comm
  real(real_8), intent(in)  :: data_in
  real(real_8), intent(inout)  :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_real8_r0'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+1*mp_double_in_bytes
  tim1=m_walltime()      
  call mpi_reduce(data_in,data_out,1,MPI_DOUBLE_PRECISION,&
       MPI_SUM,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out =data_in 
#endif
end subroutine mp_sum_root_real8_r0
subroutine mp_sum_root_real8_r1(data_in,data_out,n, root, comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root, comm
  real(real_8), intent(in) , dimension(*) :: data_in
  real(real_8), intent(inout) , dimension(*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_real8_r1'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_in_bytes
  tim1=m_walltime()      
  call mpi_reduce(data_in,data_out,n,MPI_DOUBLE_PRECISION,&
       MPI_SUM,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out (1:n)=data_in (1:n)
#endif
end subroutine mp_sum_root_real8_r1
subroutine mp_sum_root_real8_r2(data_in,data_out,n, root, comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root, comm
  real(real_8), intent(in) , dimension(1,*) :: data_in
  real(real_8), intent(inout) , dimension(1,*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_real8_r2'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_in_bytes
  tim1=m_walltime()      
  call mpi_reduce(data_in,data_out,n,MPI_DOUBLE_PRECISION,&
       MPI_SUM,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out (1,1:n)=data_in (1,1:n)
#endif
end subroutine mp_sum_root_real8_r2
subroutine mp_sum_root_complex8_r0(data_in,data_out, root, comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) ::  root, comm
  complex(real_8), intent(in)  :: data_in
  complex(real_8), intent(inout)  :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_complex8_r0'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+1*mp_double_complex_in_bytes
  tim1=m_walltime()      
  call mpi_reduce(data_in,data_out,1,MPI_DOUBLE_COMPLEX,&
       MPI_SUM,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out =data_in 
#endif
end subroutine mp_sum_root_complex8_r0
subroutine mp_sum_root_complex8_r1(data_in,data_out,n, root, comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root, comm
  complex(real_8), intent(in) , dimension(*) :: data_in
  complex(real_8), intent(inout) , dimension(*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_complex8_r1'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_complex_in_bytes
  tim1=m_walltime()      
  call mpi_reduce(data_in,data_out,n,MPI_DOUBLE_COMPLEX,&
       MPI_SUM,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out (1:n)=data_in (1:n)
#endif
end subroutine mp_sum_root_complex8_r1
subroutine mp_sum_root_complex8_r2(data_in,data_out,n, root, comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root, comm
  complex(real_8), intent(in) , dimension(1,*) :: data_in
  complex(real_8), intent(inout) , dimension(1,*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_complex8_r2'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_complex_in_bytes
  tim1=m_walltime()      
  call mpi_reduce(data_in,data_out,n,MPI_DOUBLE_COMPLEX,&
       MPI_SUM,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out (1,1:n)=data_in (1,1:n)
#endif
end subroutine mp_sum_root_complex8_r2
