MODULE rnlsm_utils
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nlm
  USE rnlsm1_utils,                    ONLY: rnlsm1
  USE rnlsm2_utils,                    ONLY: rnlsm2
  USE rnlsmd_utils,                    ONLY: rnlsmd
  USE system,                          ONLY: cntl
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlsm

CONTAINS

  ! ==================================================================
  SUBROUTINE rnlsm(c0,nstate,ikpt,ikind,tfor)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(:,:)
    INTEGER                                  :: nstate, ikpt, ikind
    LOGICAL                                  :: tfor

    CHARACTER(*), PARAMETER                  :: procedureN = 'rnlsm'

    INTEGER                                  :: isub

! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF (nlm.EQ.0) RETURN
    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    IF (cntl%tfdist) THEN
       CALL rnlsmd(c0,nstate,ikind)
    ELSE
       CALL rnlsm1(c0,nstate,ikind)
    ENDIF
    IF (tfor) CALL rnlsm2(c0,nstate,ikpt,ikind)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)

  END SUBROUTINE rnlsm
  ! ==================================================================


END MODULE rnlsm_utils
