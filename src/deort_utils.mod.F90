#include "cpmd_global.h"

MODULE deort_utils
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean,&
                                             symmat_pack,&
                                             symmat_unpack
  USE rgs_utils,                       ONLY: uinv
  USE rotate_utils,                    ONLY: rottr
  USE summat_utils,                    ONLY: summat
  USE mp_interface,                    ONLY: mp_bcast

#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: deort

CONTAINS


  ! ==================================================================
  SUBROUTINE deort(nstate,c0)
    ! ==--------------------------------------------------------------==
    ! ==    DEORTHOGONALIZE WAVEFUNCTIONS FOR VANDERBILT PP           ==
    ! ==    Transform <C0|S|C0>=1 into <C0|C0>=1                      ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(INOUT) __CONTIGUOUS &
                                             :: c0(:,:)
    INTEGER                                  :: nstate

    CHARACTER(*), PARAMETER                  :: procedureN = 'deort'

    INTEGER                                  :: ierr, isub
    INTEGER(int_8)                           :: il_smatpacked(1), il_smat(2)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: smatpacked(:), smat(:,:)
#else
    REAL(real_8), ALLOCATABLE                :: smatpacked(:), smat(:,:)
#endif
! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)
    IF(cntl%tlsd) THEN
       il_smatpacked=spin_mod%nsup*(spin_mod%nsup+1)/2+&
            spin_mod%nsdown*(spin_mod%nsdown+1)/2
    ELSE
       il_smatpacked=nstate*(nstate+1)/2
    END IF
    il_smat(1)=nstate
    il_smat(2)=nstate
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_smat,smat,procedureN//'_smat',ierr)
#else
    ALLOCATE(smat(il_smat(1),il_smat(2)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_smatpacked,smatpacked,procedureN//'_smatpacked',ierr)
#else
    ALLOCATE(smatpacked(il_smatpacked(1)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    
    CALL deort_work(nstate,c0,smat,smatpacked)

#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_smat,smat,procedureN//'_smat',ierr)
#else
    DEALLOCATE(smat,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_smatpacked,smatpacked,procedureN//'_smatpacked',ierr)
#else
    DEALLOCATE(smatpacked,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE deort
  ! ==================================================================

  SUBROUTINE deort_work(nstate,c0,smat,smatpacked)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(INOUT) __CONTIGUOUS &
                                             :: c0(:,:)
    INTEGER,INTENT(IN)                       :: nstate
    REAL(real_8),INTENT(OUT)                 :: smat(nstate,nstate)
    REAL(real_8),INTENT(OUT),CONTIGUOUS      :: smatpacked(:)

    INTEGER                                  :: i,j
    REAL(real_8)                             :: serr

! ==--------------------------------------------------------------==

    CALL ovlap(nstate,smat,c0,c0,redist=.FALSE.,full=.FALSE.)
    CALL summat(smat,nstate,symmetrization=.FALSE.,lsd=.TRUE.,gid=parai%cp_grp,&
         parent=.TRUE.)

    IF(paral%io_parent)THEN
       serr=0.0_real_8
       IF(cntl%tlsd)THEN
          !$omp parallel private(i,j)reduction(+:serr)
          !$omp do
          DO i=1,spin_mod%nsup
             DO j=1,i
                serr=serr+smat(j,i)
             ENDDO
          ENDDO
          !$omp end do nowait
          !$omp do
          DO i=spin_mod%nsup+1,nstate
             DO j=spin_mod%nsup+1,i
                serr=serr+smat(j,i)
             ENDDO
          ENDDO
          !$omp end parallel
       ELSE
          !$omp parallel do private(i,j) reduction(+:serr)
          DO i=1,nstate
             DO j=1,i
                serr=serr+smat(j,i)
             END DO
          END DO
       END IF
       serr=serr-REAL(nstate,kind=real_8)
    END IF
    CALL mp_bcast(serr,parai%io_source,parai%cp_grp)
    IF (ABS(serr).GT.1.e-8_real_8) THEN    
       IF(paral%io_parent)THEN
          IF (cntl%tlsd) THEN
             CALL uinv('U',smat(1,1),nstate,spin_mod%nsup)
             CALL uinv('U',smat(spin_mod%nsup+1,spin_mod%nsup+1),nstate,spin_mod%nsdown)
             CALL symmat_pack(smat,smatpacked,nstate,spin_mod%nsup,spin_mod%nsdown)
          ELSE
             CALL uinv('U',smat,nstate,nstate)
             CALL symmat_pack(smat,smatpacked,nstate,nstate,0)
          ENDIF
       END IF
       CALL mp_bcast(smatpacked,SIZE(smatpacked),parai%io_source,parai%cp_grp)
       IF(cntl%tlsd)THEN
          CALL symmat_unpack(smat,smatpacked,nstate,spin_mod%nsup,spin_mod%nsdown,.FALSE.)
       ELSE
          CALL symmat_unpack(smat,smatpacked,nstate,nstate,0,.FALSE.)
       END IF
       CALL rottr(1._real_8,c0,smat,"N",nstate,ncpw%ngw,cntl%tlsd,spin_mod%nsup,&
            spin_mod%nsdown,use_cp=.TRUE.,redist=.TRUE.)
    END IF
    IF (geq0) CALL zclean(c0,nstate,ncpw%ngw)
    RETURN
  END SUBROUTINE deort_work
  ! ==================================================================

END MODULE deort_utils
