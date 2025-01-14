#include "cpmd_global.h"

MODULE cofor_utils
  USE cppt,                            ONLY: gk,&
                                             inyh
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes
#ifdef _VERBOSE_FORCE_DBG
  USE error_handling,                  ONLY: stopgm
#endif
  USE ions,                            ONLY: ions1,&
                                             ions0
  USE mp_interface,                    ONLY: mp_sum
  USE kinds,                           ONLY: real_8
  USE nlcc,                            ONLY: corel,&
                                             rhoc,&
                                             vnlcc,&
                                             vnlt
#ifdef _VERBOSE_FORCE_DBG
  USE parac,                           ONLY: parai, &
                                             paral
#else
  USE parac,                           ONLY: parai
#endif
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             ncpw,&
                                             parm,&
                                             maxsys
  USE timer,                           ONLY: tihalt,&
                                             tiset


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cofor

CONTAINS

  ! ==================================================================
  SUBROUTINE cofor(fion,vpot)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == FORCES ON IONS DUE TO CORE CHARGES                           ==
    ! ==--------------------------------------------------------------==
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: fion(:,:,:)
    COMPLEX(real_8)                          :: vpot(ncpw%nhg,*)

    COMPLEX(real_8)                          :: ei123, vxc
    INTEGER                                  :: ia, ig, is, isa, k, isa0, &
                                                ig_start, ig_end, isub
    REAL(real_8)                             :: omtp, vcgs, fiont(3)
#ifdef _VERBOSE_FORCE_DBG
    INTEGER                                  :: ierr
    REAL(real_8),ALLOCATABLE                 :: dbg_forces(:,:,:)
#endif
    CHARACTER(*), PARAMETER                  :: procedureN = 'cofor'

    CALL tiset(procedureN,isub)

    IF(parai%cp_nogrp.GT.1.AND.parai%cp_inter_me.GT.0)THEN
       !$omp parallel do private(is,ia)
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             fion(1:3,ia,is)=0._real_8
          END DO
       END DO
    END IF
    CALL cp_grp_get_sizes(first_nhg=ig_start,last_nhg=ig_end)
    omtp=2._real_8*parm%omega*parm%tpiba
    IF(cntl%bigmem)THEN
       IF(cntl%tlsd)THEN
          !$omp parallel private(isa0,is,ia,fiont,isa,ig,vxc,vcgs,k)
          isa0=0
          DO is=1,ions1%nsp
             IF(corel%tnlcc(is))THEN
                !$omp do
                DO ia=1,ions0%na(is)
                   fiont=0._real_8
                   isa=isa0+ia
                   DO ig=ig_start,ig_end
                      vxc=0.5_real_8*(vpot(ig,1)+vpot(ig,2)-2*vnlt(ig)-&
                           vnlcc(ig,1)-vnlcc(ig,2))
                      vcgs=-AIMAG(CONJG(vxc)*eigrb(ig,isa)*rhoc(ig,is))
                      DO k=1,3
                         fiont(k)=fiont(k)+gk(k,ig)*vcgs*omtp
                      END DO
                   END DO
                   DO k=1,3
                      fion(k,ia,is)=fion(k,ia,is)+fiont(k)
                   END DO
                END DO
                !$omp end do nowait
             END IF
             isa0=isa0+ions0%na(is)
          END DO
          !$omp end parallel
       ELSE
          !$omp parallel private(isa0,is,ia,fiont,isa,ig,vxc,vcgs,k)
          isa0=0
          DO is=1,ions1%nsp
             IF(corel%tnlcc(is))THEN
                !$omp do
                DO ia=1,ions0%na(is)
                   fiont=0._real_8
                   isa=isa0+ia
                   DO ig=ig_start,ig_end
                      vxc=vpot(ig,1)-vnlt(ig)-vnlcc(ig,1)
                      vcgs=-AIMAG(CONJG(vxc)*eigrb(ig,isa)*rhoc(ig,is))
                      DO k=1,3
                         fiont(k)=fiont(k)+gk(k,ig)*vcgs*omtp
                      END DO
                   END DO
                   DO k=1,3
                      fion(k,ia,is)=fion(k,ia,is)+fiont(k)
                   END DO
                END DO
                !$omp end do nowait
             END IF
             isa0=isa0+ions0%na(is)
          END DO
          !$omp end parallel
       END IF
    ELSE
       !$omp parallel do private(isa,ia,is,ig,ei123,vxc,vcgs)
       DO isa=1,ions1%nat
          ia=iatpt(1,isa)
          is=iatpt(2,isa)
          IF (corel%tnlcc(is)) THEN
             DO ig=ig_start,ig_end
                ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))
                vxc=vpot(ig,1)-vnlt(ig)-vnlcc(ig,1)
                IF (cntl%tlsd) vxc=0.5_real_8*(vxc+vpot(ig,2)-vnlt(ig)-vnlcc(ig,2))
                vcgs=-AIMAG(CONJG(vxc)*ei123*rhoc(ig,is))
                fion(1,ia,is)=fion(1,ia,is)+gk(1,ig)*vcgs*omtp
                fion(2,ia,is)=fion(2,ia,is)+gk(2,ig)*vcgs*omtp
                fion(3,ia,is)=fion(3,ia,is)+gk(3,ig)*vcgs*omtp
             ENDDO
          ENDIF
       ENDDO
    END IF

    IF(parai%cp_nogrp.GT.1)THEN
       CALL mp_sum(fion,3*maxsys%nax*maxsys%nsx,parai%cp_inter_grp)
    END IF

#ifdef _VERBOSE_FORCE_DBG
    ALLOCATE(dbg_forces(3,maxsys%nax,maxsys%nsx), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate dbg_forces',& 
         __LINE__,__FILE__)
    dbg_forces=fion
    CALL mp_sum(dbg_forces,3*maxsys%nax*maxsys%nsx,parai%allgrp)
    IF (paral%io_parent) THEN
       WRITE(6,*) "===================================="
       WRITE(6,*) "DEBUG FORCES", procedureN
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             WRITE(6,*) dbg_forces(1:3,ia,is),ia,is
          END DO
       END DO
    END IF
    DEALLOCATE(dbg_forces,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#endif
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cofor
  ! ==================================================================

END MODULE cofor_utils
