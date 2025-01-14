#include "cpmd_global.h"

MODULE autotune_utils
  USE elct,                            ONLY: crge
  USE spin,                            ONLY: clsd
  USE rhoofr_utils,                    ONLY: rhoofr_batchfft
  USE vpsi_utils,                      ONLY: vpsi_batchfft
  USE rnlsm_utils,                     ONLY: rnlsm
  USE fft,                             ONLY: batch_fft,&
                                             fft_tune_max_it
  USE fftprp_utils,                    ONLY: autotune_fftbatchsize
  USE system,                          ONLY: cnti, cntl
  USE rswfmod,                         ONLY: rsactive
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC autotune
  
CONTAINS
  ! ==================================================================
  SUBROUTINE autotune(c0,c2,rhoe,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! == Run all iterations for autotuning rnlsm/vpsi/rhoofr          ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(IN) __CONTIGUOUS  :: c0(:,:)
    REAL(real_8),INTENT(OUT) __CONTIGUOUS    :: rhoe(:,:)
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: psi(:,:), c2(:,:)
    INTEGER, INTENT(IN)                      :: nstate

    CHARACTER(*), PARAMETER                  :: procedureN = 'autotuning'

    INTEGER                                  :: isub, it, num_it

    CALL tiset(procedureN,isub)
    
    IF(cntl%fft_tune_batchsize.OR.cntl%rnlsm_autotune)THEN
       num_it=MAX(fft_tune_max_it,cnti%rnlsm_autotune_maxit)
       DO it=1,num_it
          IF(it.LE.cnti%rnlsm_autotune_maxit.AND.cntl%overlapp_comm_comp)THEN
             CALL rnlsm(c0,nstate,1,1,.FALSE.,unpack_dfnl_fnl=.FALSE.)
          END IF
          IF(it.LE.fft_tune_max_it.AND.batch_fft)THEN
             rsactive = cntl%krwfn
             CALL autotune_fftbatchsize()
             CALL rhoofr_batchfft(c0,rhoe,psi(:,1),nstate)
             CALL vpsi_batchfft(c0,c2,crge%f(:,1),rhoe,psi(:,1),nstate,1,clsd%nlsd,.TRUE.)
             rsactive = .FALSE.
          END IF
       END DO
    END IF
    CALL autotune_fftbatchsize()
    cntl%rnlsm_autotune=.FALSE.
    CALL tihalt(procedureN,isub)
  END SUBROUTINE autotune

END MODULE autotune_utils
