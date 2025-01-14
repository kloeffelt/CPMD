#include "cpmd_global.h"
MODULE potfor_utils
  USE cppt,                            ONLY: gk,&
                                             inyh,&
                                             nzh,&
                                             rhops,&
                                             scg,&
                                             vps
#ifdef _VERBOSE_FORCE_DBG
  USE error_handling,                  ONLY: stopgm
#endif
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
#ifdef _VERBOSE_FORCE_DBG 
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
#endif
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
#ifdef _VERBOSE_FORCE_DBG
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             iatpt,&
                                             ncpw,&
                                             parm,&
                                             maxsys
#else
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             iatpt,&
                                             ncpw,&
                                             parm
#endif
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: potfor
  PUBLIC :: potabfor

CONTAINS

  ! ==================================================================
  SUBROUTINE potfor(fion,v,eirop)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == POTENTIAL ENERGY CONTRIBUTIONS TO THE FORCES ON THE IONS     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)
    COMPLEX(real_8)                          :: v(*), eirop(ncpw%nhg)

    CHARACTER(*), PARAMETER                  :: procedureN = 'potfor'

    COMPLEX(real_8)                          :: ei123, gx, gy, gz, rhet, &
                                                rhets, rhog, rhogs, rp, txx, &
                                                tyy, tzz, vcgs
    INTEGER                                  :: ia, ig, ig1, is, isa, isub
    REAL(real_8)                             :: omtp, ft1, ft2, ft3
#ifdef _VERBOSE_FORCE_DBG
    INTEGER                                  :: ierr
    REAL(real_8),ALLOCATABLE                 :: dbg_forces(:,:,:)
#endif

    CALL tiset(procedureN,isub)
    omtp=2._real_8*parm%omega*parm%tpiba
    ig1=1
    IF (geq0) ig1=2
    IF (cntl%bigmem) THEN
       !$omp parallel do private(ISA,IA,IS,IG,EI123,RP,RHET,RHOG,RHETS,RHOGS, &
       !$omp  GX,GY,GZ,VCGS,ft1,ft2,ft3)
       DO isa=1,ions1%nat
          ia=iatpt(1,isa)
          is=iatpt(2,isa)
          ft1=0._real_8
          ft2=0._real_8
          ft3=0._real_8
          DO ig=ig1,ncpw%nhg
             ei123=eigrb(ig,isa)
             rp=eirop(ig)
             rhet=v(nzh(ig))
             rhog=rhet+rp
             rhets=CONJG(rhet)
             rhogs=CONJG(rhog)
             gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
             gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
             gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
             vcgs=scg(ig)*rhogs
             ei123=ei123*(rhops(is,ig)*vcgs+vps(is,ig)*rhets)
             ft1=ft1+REAL(ei123*gx)
             ft2=ft2+REAL(ei123*gy)
             ft3=ft3+REAL(ei123*gz)
          ENDDO
          fion(1,ia,is)=fion(1,ia,is)+ft1*omtp
          fion(2,ia,is)=fion(2,ia,is)+ft2*omtp
          fion(3,ia,is)=fion(3,ia,is)+ft3*omtp
       ENDDO
    ELSE
       !$omp parallel do private(ISA,IA,IS,IG,EI123,RP,RHET,RHOG,RHETS,RHOGS, &
       !$omp  GX,GY,GZ,VCGS,ft1,ft2,ft3)
       DO isa=1,ions1%nat
          ia=iatpt(1,isa)
          is=iatpt(2,isa)
          DO ig=ig1,ncpw%nhg
             ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                  ei3(isa,inyh(3,ig))
             rp=eirop(ig)
             rhet=v(nzh(ig))
             rhog=rhet+rp
             rhets=CONJG(rhet)
             rhogs=CONJG(rhog)
             gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
             gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
             gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
             vcgs=scg(ig)*rhogs
             ei123=ei123*(rhops(is,ig)*vcgs+vps(is,ig)*rhets)
             ft1=ft1+REAL(ei123*gx)
             ft2=ft2+REAL(ei123*gy)
             ft3=ft3+REAL(ei123*gz)
          ENDDO
          fion(1,ia,is)=fion(1,ia,is)+ft1*omtp
          fion(2,ia,is)=fion(2,ia,is)+ft2*omtp
          fion(3,ia,is)=fion(3,ia,is)+ft3*omtp
       ENDDO
    ENDIF
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
  END SUBROUTINE potfor
  ! ==================================================================
  SUBROUTINE potabfor(fion,rhog)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == POTENTIAL ENERGY CONTRIBUTIONS TO THE FORCES ON THE IONS     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)
    COMPLEX(real_8)                          :: rhog(fpar%nnr1)

    CHARACTER(*), PARAMETER                  :: procedureN = 'potabfor'

    COMPLEX(real_8)                          :: ei123, gx, gy, gz, rg, txx, &
                                                tyy, tzz
    INTEGER                                  :: ia, ig, ig1, is, isa, isub
    REAL(real_8)                             :: omtp

    CALL tiset(procedureN,isub)
    omtp=2._real_8*parm%omega*parm%tpiba
    ig1=1
    IF (geq0) ig1=2
    isa=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          isa=isa+1
          IF (cntl%bigmem) THEN
             !$omp parallel do private(IG,EI123,RG,GX,GY,GZ,TXX,TYY,TZZ) &
             !$omp  shared(ISA,RHOG)
             DO ig=ig1,ncpw%nhg
                ei123=eigrb(ig,isa)
                rg=CONJG(rhog(ig))
                ! mb              GX=CMPLX(0._real_8,GK(1,IG))
                ! mb              GY=CMPLX(0._real_8,GK(2,IG))
                ! mb              GZ=CMPLX(0._real_8,GK(3,IG))
                txx=vps(is,ig)*rg*CMPLX(0._real_8,gk(1,ig),kind=real_8)
                tyy=vps(is,ig)*rg*CMPLX(0._real_8,gk(2,ig),kind=real_8)
                tzz=vps(is,ig)*rg*CMPLX(0._real_8,gk(3,ig),kind=real_8)
                fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*txx)*omtp
                fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*tyy)*omtp
                fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*tzz)*omtp
             ENDDO
          ELSE
             !$omp parallel do private(IG,EI123,RG,GX,GY,GZ,TXX,TYY,TZZ) &
             !$omp  shared(ISA,RHOG)
             DO ig=ig1,ncpw%nhg
                ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))
                rg=CONJG(rhog(ig))
                ! mb              GX=CMPLX(0._real_8,GK(1,IG))
                ! mb              GY=CMPLX(0._real_8,GK(2,IG))
                ! mb              GZ=CMPLX(0._real_8,GK(3,IG))
                txx=vps(is,ig)*rg*CMPLX(0._real_8,gk(1,ig),kind=real_8)
                tyy=vps(is,ig)*rg*CMPLX(0._real_8,gk(2,ig),kind=real_8)
                tzz=vps(is,ig)*rg*CMPLX(0._real_8,gk(3,ig),kind=real_8)
                fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*txx)*omtp
                fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*tyy)*omtp
                fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*tzz)*omtp
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE potabfor
  ! ==================================================================

END MODULE potfor_utils
