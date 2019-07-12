module scratch_interface
  use scratch_data, only : dp, int_long
  use, intrinsic :: iso_c_binding
  use scratch_segment_managment_interface
  implicit none
  private
  public :: request_scratch
  public :: save_scratch
  public :: free_scratch
  public :: request_saved_scratch

  interface request_scratch
     module procedure :: request_scratch_r1
     module procedure :: request_scratch_r2
     module procedure :: request_scratch_r3
     module procedure :: request_scratch_r4
     module procedure :: request_scratch_r5
     module procedure :: request_scratch_r6
     module procedure :: request_scratch_c1
     module procedure :: request_scratch_c2
     module procedure :: request_scratch_c3
     module procedure :: request_scratch_c4
     module procedure :: request_scratch_c5
     module procedure :: request_scratch_c6
  end interface request_scratch

  interface save_scratch
     module procedure :: save_scratch_r1
     module procedure :: save_scratch_r2
     module procedure :: save_scratch_r3
     module procedure :: save_scratch_r4
     module procedure :: save_scratch_r5
     module procedure :: save_scratch_r6
     module procedure :: save_scratch_c1
     module procedure :: save_scratch_c2
     module procedure :: save_scratch_c3
     module procedure :: save_scratch_c4
     module procedure :: save_scratch_c5
     module procedure :: save_scratch_c6
  end interface save_scratch

  interface free_scratch
     module procedure :: free_scratch_r1
     module procedure :: free_scratch_r2
     module procedure :: free_scratch_r3
     module procedure :: free_scratch_r4
     module procedure :: free_scratch_r5
     module procedure :: free_scratch_r6
     module procedure :: free_scratch_c1
     module procedure :: free_scratch_c2
     module procedure :: free_scratch_c3
     module procedure :: free_scratch_c4
     module procedure :: free_scratch_c5
     module procedure :: free_scratch_c6
  end interface free_scratch

  interface request_saved_scratch
     module procedure :: request_saved_scratch_r1
     module procedure :: request_saved_scratch_r2
     module procedure :: request_saved_scratch_r3
     module procedure :: request_saved_scratch_r4
     module procedure :: request_saved_scratch_r5
     module procedure :: request_saved_scratch_r6
     module procedure :: request_saved_scratch_c1
     module procedure :: request_saved_scratch_c2
     module procedure :: request_saved_scratch_c3
     module procedure :: request_saved_scratch_c4
     module procedure :: request_saved_scratch_c5
     module procedure :: request_saved_scratch_c6
  end interface request_saved_scratch

  contains
include "scratch_interface_templates.f90"

end module scratch_interface
