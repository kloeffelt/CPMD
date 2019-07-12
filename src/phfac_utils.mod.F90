#include "cpmd_global.h"

MODULE phfac_utils
  USE cppt,                            ONLY: inyh
  USE error_handling,                  ONLY: stopgm
  USE gvec,                            ONLY: gvec_com
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE kpnt,                            ONLY: eigkr,&
                                             eikr,&
                                             rk
  USE kpts,                            ONLY: tkpts
  !$ USE omp_lib,                         ONLY: omp_get_thread_num
  USE parac,                           ONLY: paral, &
                                             parai
  USE prmem_utils,                     ONLY: prmem
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigr,&
                                             eigrb,&
                                             natx
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             kpbeg,&
                                             ncpw,&
                                             nkpbl,&
                                             nkpt,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: phfac
  PUBLIC :: calc_eigkr

CONTAINS

  ! ==================================================================
  SUBROUTINE phfac(tau0)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'phfac'

    COMPLEX(real_8)                          :: ctem1, ctem2, ctem3, ctep1, &
                                                ctep2, ctep3, ei10, ei20, &
                                                ei30, svtmpm, svtmpp, zsum
    INTEGER                                  :: i, ia, ierr, ig, ik, ikk, &
                                                ikpt, is, isa, isub, j, k, &
                                                nh1, nh2, nh3,  methread
    INTEGER(int_8)                           :: il_ei2t(2), il_ei3t(2), il_ei1t(2)
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: ar1, ar2, ar3, sum, sum1, &
                                                sum2, sum3
#ifdef _USE_SCRATCHLIBRARY
    COMPLEX(real_8),POINTER __CONTIGUOUS     :: ei1t(:,:), ei2t(:,:), ei3t(:,:)
#else
    COMPLEX(real_8),ALLOCATABLE              :: ei1t(:,:), ei2t(:,:), ei3t(:,:)
#endif
! ==--------------------------------------------------------------==

    IF (spar%nr1s.LT.3) THEN
       CALL stopgm('PHFAC',' PHFAC: NR1 TOO SMALL ',& 
            __LINE__,__FILE__)
    ENDIF
    IF (spar%nr2s.LT.3) THEN
       CALL stopgm('PHFAC',' PHFAC: NR2 TOO SMALL ',& 
            __LINE__,__FILE__)
    ENDIF
    IF (spar%nr3s.LT.3) THEN
       CALL stopgm('PHFAC',' PHFAC: NR3 TOO SMALL ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)
    IF (ifirst.EQ.0) THEN
       !TK either eigrb or ei1-3
       IF (cntl%bigmem) THEN
          ALLOCATE(eigrb(ncpw%nhg,ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(ei1(natx,(2*spar%nr1s-1)),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(ei2(natx,(2*spar%nr2s-1)),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(ei3(natx,(2*spar%nr3s-1)),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ALLOCATE(eigr(ncpw%ngw,ions1%nat,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)! FIXME deallocate missing
       IF (tkpts%tkpnt) THEN
          ALLOCATE(eikr(nkpt%nkpts,ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)

          ALLOCATE(eigkr(nkpt%ngwk,ions1%nat,nkpt%nkpnt),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ifirst  = 1
       IF (paral%parent) THEN
          CALL prmem('     PHFAC')
       ENDIF
    ENDIF
    !TK temporary arrays
    il_ei1t(1)=2*spar%nr1s-1
    il_ei1t(2)=parai%ncpus
    il_ei2t(1)=2*spar%nr2s-1
    il_ei2t(2)=parai%ncpus
    il_ei3t(1)=2*spar%nr3s-1
    il_ei3t(2)=parai%ncpus
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_ei1t,ei1t,procedureN//'_ei1t')
    CALL request_scratch(il_ei2t,ei2t,procedureN//'_ei2t')
    CALL request_scratch(il_ei3t,ei3t,procedureN//'_ei3t')
#else
    ALLOCATE(ei1t(il_ei1t(1),il_ei2t(2)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ei2t(il_ei2t(1),il_ei2t(2)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ei3t(il_ei3t(1),il_ei3t(2)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
#endif
    ! ==--------------------------------------------------------------==
    nh1=spar%nr1s/2
    nh2=spar%nr2s/2
    nh3=spar%nr3s/2
    methread=1
    !$omp parallel private(isa,ia,is,sum1,sum2,sum3,ar1,ar2,ar3) &
    !$omp private(ctep1,ctep2,ctep3,ctem1,ctem2,ctem3) &
    !$omp private(svtmpp,svtmpm,i,j,k,ei10,ei20,ei30,methread) &
    !$omp shared(nh1,nh2,nh3)
    !$ methread=omp_get_thread_num()+1
    !$omp do
    DO isa=1,ions1%nat
       ia=iatpt(1,isa)
       is=iatpt(2,isa)
       sum1=gvec_com%b1(1)*tau0(1,ia,is)+gvec_com%b1(2)*tau0(2,ia,is)+&
            gvec_com%b1(3)*tau0(3,ia,is)
       sum2=gvec_com%b2(1)*tau0(1,ia,is)+gvec_com%b2(2)*tau0(2,ia,is)+&
            gvec_com%b2(3)*tau0(3,ia,is)
       sum3=gvec_com%b3(1)*tau0(1,ia,is)+gvec_com%b3(2)*tau0(2,ia,is)+&
            gvec_com%b3(3)*tau0(3,ia,is)
       ar1=parm%tpiba*sum1
       ar2=parm%tpiba*sum2
       ar3=parm%tpiba*sum3
       ctep1=CMPLX(COS(ar1),-SIN(ar1),kind=real_8)
       ctep2=CMPLX(COS(ar2),-SIN(ar2),kind=real_8)
       ctep3=CMPLX(COS(ar3),-SIN(ar3),kind=real_8)
       ctem1=CONJG(ctep1)
       ctem2=CONJG(ctep2)
       ctem3=CONJG(ctep3)

       ei10=ctep1**(-nh1)
       ei1t(1,methread)=CMPLX(1.0_real_8,0.0_real_8,kind=real_8)*ei10
       svtmpp=ctep1
       svtmpm=ctem1
       DO i=2,spar%nr1s
          ei1t(i,methread)=svtmpp*ei10
          svtmpp=svtmpp*ctep1
          ei1t(spar%nr1s+i-1,methread)=svtmpm*ei10
          svtmpm=svtmpm*ctem1
       END DO

       ei20=ctep2**(-nh2)
       ei2t(1,methread)=CMPLX(1.0_real_8,0.0_real_8,kind=real_8)*ei20
       svtmpp=ctep2
       svtmpm=ctem2
       DO j=2,spar%nr2s
          ei2t(j,methread)=svtmpp*ei20
          svtmpp=svtmpp*ctep2
          ei2t(spar%nr2s+j-1,methread)=svtmpm*ei20
          svtmpm=svtmpm*ctem2
       END DO

       ei30=ctep3**(-nh3)
       ei3t(1,methread)=CMPLX(1.0_real_8,0.0_real_8,kind=real_8)*ei30
       svtmpp=ctep3
       svtmpm=ctem3
       DO k=2,spar%nr3s
          ei3t(k,methread)=svtmpp*ei30
          svtmpp=svtmpp*ctep3
          ei3t(spar%nr3s+k-1,methread)=svtmpm*ei30
          svtmpm=svtmpm*ctem3
       END DO

       ! ==--------------------------------------------------------------==
       IF (cntl%bigmem) THEN
          DO ig=1,ncpw%ngw
             eigr(ig,isa,1)=ei1t(inyh(1,ig),methread)*ei2t(inyh(2,ig),methread)&
                  *ei3t(inyh(3,ig),methread)
             eigrb(ig,isa)=ei1t(inyh(1,ig),methread)*ei2t(inyh(2,ig),methread)&
                  *ei3t(inyh(3,ig),methread)
          END DO
          DO ig=ncpw%ngw+1,ncpw%nhg
             eigrb(ig,isa)=ei1t(inyh(1,ig),methread)*ei2t(inyh(2,ig),methread)&
                  *ei3t(inyh(3,ig),methread)
          END DO
       ELSE
       !TK ei1-3 only needed if bigmem is not active
          DO ig=1,ncpw%ngw
             eigr(ig,isa,1)=ei1t(inyh(1,ig),methread)*ei2t(inyh(2,ig),methread)&
                  *ei3t(inyh(3,ig),methread)
          END DO
          DO ig=1,2*spar%nr1s-1
             ei1(isa,ig)=ei1t(ig,methread)
          END DO
          DO ig=1,2*spar%nr2s-1
             ei2(isa,ig)=ei2t(ig,methread)
          END DO
          DO ig=1,2*spar%nr3s-1
             ei3(isa,ig)=ei3t(ig,methread)
          END DO
       END IF
    END DO
    !$omp end do nowait
    !$omp end parallel
    ! ==--------------------------------------------------------------==
    IF (tkpts%tkpnt) THEN
       DO ikpt=1,nkpt%nblkp
1         DO ik=1,nkpbl(ikpt)
             ikk=kpbeg(ikpt)+ik
             !$omp parallel do private(ISA,IA,IS,SUM1,SUM2,SUM3,SUM,ZSUM,IG)
             DO isa=1,ions1%nat
                ia=iatpt(1,isa)
                is=iatpt(2,isa)
                sum1=gvec_com%b1(1)*tau0(1,ia,is)+gvec_com%b1(2)*tau0(2,ia,is)+&
                     gvec_com%b1(3)*tau0(3,ia,is)
                sum2=gvec_com%b2(1)*tau0(1,ia,is)+gvec_com%b2(2)*tau0(2,ia,is)+&
                     gvec_com%b2(3)*tau0(3,ia,is)
                sum3=gvec_com%b3(1)*tau0(1,ia,is)+gvec_com%b3(2)*tau0(2,ia,is)+&
                     gvec_com%b3(3)*tau0(3,ia,is)
                sum=rk(1,ikk)*sum1+rk(2,ikk)*sum2+rk(3,ikk)*sum3
                zsum=CMPLX(COS(sum),SIN(sum),kind=real_8)
                eikr(ikk,isa)=zsum
                DO ig=1,ncpw%ngw
                   eigkr(ig,isa,ik)=eigr(ig,isa,1)*zsum
                   eigkr(ig+ncpw%ngw,isa,ik)=CONJG(eigr(ig,isa,1))*zsum
                END DO
             END DO
          END DO
          IF (tkpts%tkblock) THEN
             CALL wkpt_swap(eigkr,1,ikpt,'EIGKR')
          ENDIF
       END DO
    ELSE
       ! kpt..not to duplicate code.
       ! NB: it is impossible to achieve this using Fortran 90, because
       ! EIGKR has rank 3, and EIGR has rank 2. Even if 
       ! EIGKR has dimensions (:,:,1) it cannot point to (:,:)
       ! To refactor this we will have to either make EIGR rank 3 with last dim = 1,
       ! or allocate memory for EIGKR and copy...
       ! 
       ! semantically this should be the following
       ! EIGKR => EIGR
       ! in principle, EIGKR should be just removed (?), and EIGR should have rank 3
       ! where last dim is NKPNT

       ! TODO: fix this!
       eigkr => eigr
    ENDIF

#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_ei3t,ei3t,procedureN//'_ei3t')
    CALL free_scratch(il_ei2t,ei2t,procedureN//'_ei2t')
    CALL free_scratch(il_ei1t,ei1t,procedureN//'_ei1t')
#else
    DEALLOCATE(ei1t,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ei2t,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ei3t,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
#endif

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE phfac
  ! ==================================================================
  SUBROUTINE calc_eigkr(ikpt)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES EIGKR (USED WITH THE OPTION TKBLOCK)              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ikpt

    COMPLEX(real_8)                          :: zsum
    INTEGER                                  :: ia, ig, ik, ikind, is, isa, &
                                                isub

    CALL tiset('CALC_EIGKR',isub)
    DO ikind=1,nkpbl(ikpt)
       ik=kpbeg(ikpt)+ikind
       isa=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             isa=isa+1
             zsum=eikr(ik,isa)
             DO ig=1,ncpw%ngw
                eigkr(ig,isa,ikind)=eigr(ig,isa,1)*zsum
                eigkr(ig+ncpw%ngw,isa,ikind)=CONJG(eigr(ig,isa,1))*zsum
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CALL tihalt('CALC_EIGKR',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE calc_eigkr
  ! ==================================================================

END MODULE phfac_utils
