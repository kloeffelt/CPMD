MODULE rekine_utils
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes
  USE dotp_utils,                      ONLY: dotp_c1_cp
  USE geq0mod,                         ONLY: geq0
  USE harm,                            ONLY: xmu
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rekine

CONTAINS

  ! ==================================================================
  SUBROUTINE rekine(cm,nstate,ekinc,use_cp_grps)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate
    COMPLEX(real_8),INTENT(IN)               :: cm(ncpw%ngw,nstate)
    REAL(real_8),INTENT(OUT)                 :: ekinc
    LOGICAL,INTENT(IN),OPTIONAL              :: use_cp_grps

    INTEGER                                  :: i, ig, isub, ngw_local,&
                                                ibeg_c0, iend_c0, gid
    LOGICAL                                  :: geq0_local, cp_active
    REAL(real_8)                             :: ax, bx, pf
    CHARACTER(*), PARAMETER                  :: procedureN = 'rekine'
! ==--------------------------------------------------------------==
! ==  COMPUTE FICTITIOUS KINETIC ENERGY OF THE ELECTRONS          ==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF(PRESENT(use_cp_grps))THEN
       cp_active=use_cp_grps
    ELSE
       cp_active=.FALSE.
    END IF
    IF(cp_active)THEN
       CALL cp_grp_get_sizes(ngw_l=ngw_local,first_g=ibeg_c0,geq0_l=geq0_local)
       gid=parai%cp_grp
    ELSE
       ngw_local=ncpw%ngw
       ibeg_c0=1
       gid=parai%allgrp
       geq0_local=geq0
    END IF

    ekinc=0._real_8
    IF (cntl%tmass) THEN
       !$omp parallel do private(i,ig,pf,ax,bx) reduction(+:ekinc)
       DO i=1,nstate
          DO ig=ibeg_c0,iend_c0
             pf=2.0_real_8*xmu(ig)
             ax=REAL(cm(ig,i))
             bx=AIMAG(cm(ig,i))
             ekinc=ekinc+pf*(ax*ax+bx*bx)
          ENDDO
          IF (geq0_local) ekinc=ekinc-xmu(1)*REAL(cm(1,i))*REAL(cm(1,i))
       ENDDO
    ELSE
       !$omp parallel do private(i) reduction(+:ekinc)
       DO i=1,nstate
          ekinc=ekinc+dotp_c1_cp(ngw_local,cm(ibeg_c0,i),geq0_local)
       ENDDO
       ekinc=ekinc*cntr%emass
    ENDIF
    CALL mp_sum(ekinc,gid)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rekine
  ! ==================================================================

END MODULE rekine_utils
