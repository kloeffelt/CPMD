MODULE cvan
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == VANDERBILD PSEUDOPOTENTIALS                                  ==
  ! ==--------------------------------------------------------------==
  INTEGER, ALLOCATABLE :: nelev(:)

  REAL(real_8), ALLOCATABLE :: qq(:,:,:)
  REAL(real_8), ALLOCATABLE :: dvan(:,:,:)
  REAL(real_8), ALLOCATABLE :: deeq(:,:,:,:)
  COMPLEX(real_8), ALLOCATABLE :: qg(:,:)
  ! ==================================================================

END MODULE cvan
