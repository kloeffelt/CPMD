#include "cpmd_global.h"

MODULE newd_utils
  USE cppt,                            ONLY: gk,&
                                             inyh
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes,&
                                             cp_grp_split_atoms
  USE cvan,                            ONLY: qg
  USE distribution_utils,              ONLY: dist_entity
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: nlps_com
  !$ USE omp_lib,                       ONLY: omp_get_thread_num
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE qvan2_utils,                     ONLY: qvan2
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb,&
                                             fnl_packed,&
                                             il_fnl_packed
  USE system,                          ONLY: cntl,&
                                             iatpe,&
                                             ipept,&
                                             maxsys,&
                                             ncpw,&
                                             parm,&
                                             cnti
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: newd

CONTAINS

    ! ==--------------------------------------------------------------==
    ! ==         THIS ROUTINE CALCULATES THE ARRAY DEEQ               ==
    ! ==                                                              ==
    ! ==  DEEQ_I,IJ = OMEGA( V(G=0) Q_I,IJ(G=0) +                     ==
    ! ==       2 SUM_G> RE[ V*(G) Q_I,IJ(G) E^-IG.R_I ] )             ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
  SUBROUTINE newd(i_start,i_end,deeq,f,vpot,fion,tfor)
    REAL(REAL_8),INTENT(OUT)                 :: deeq(ions1%nat,maxsys%nhxs,*)
    COMPLEX(real_8),INTENT(IN)               :: vpot(ncpw%nhg)
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: fion(:,:,:)
    INTEGER,INTENT(IN)                       :: i_start, i_end
    REAL(real_8),INTENT(IN)                  :: f(*)
    LOGICAL,INTENT(IN)                       :: tfor
    INTEGER                                  :: is, ia, isub, nst(2,0:parai%nproc-1), &
                                                na_grp(2,ions1%nsp,0:parai%cp_nogrp-1), &
                                                na(2,ions1%nsp), ig_start, nhg_loc

    CHARACTER(*), PARAMETER                  :: procedureN='newd'

    CALL tiset(procedureN,isub)
    IF (cntl%tfdist) CALL stopgm(procedureN,'TFDIST NOT IMPLEMENTED',&
         __LINE__,__FILE__)
    CALL cp_grp_split_atoms(na_grp)
    na(:,:)=na_grp(:,:,parai%cp_inter_me)
    CALL cp_grp_get_sizes(first_nhg=ig_start,nhg_l=nhg_loc)
    !generate 'new' nst12 mapping -> nwa12
    CALL dist_entity(i_end-i_start+1,parai%nproc,nst)
    !shift to starting state
    nst(1,parai%me)=nst(1,parai%me)+i_start-1
    nst(2,parai%me)=nst(2,parai%me)+i_start-1

    IF(tfor)THEN
       IF (parai%cp_nogrp.GT.1 .AND. parai%cp_inter_me .GT. 0) THEN
          !$omp parallel do private(is,ia)
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                fion(1:3,ia,is)=0._real_8
             END DO
          END DO
       END IF
    END IF

    IF (cntl%bigmem) THEN
       CALL prep_bigmem_newd(nst(1,parai%me),ig_start,nhg_loc,na,deeq,f,fion,vpot,&
            fnl_packed,tfor)
    ELSE
       CALL prep_smallmem_newd(nst(1,parai%me),ig_start,nhg_loc,na,deeq,f,fion,vpot,&
            fnl_packed,tfor)
    ENDIF

    IF(tfor)THEN
       IF (parai%cp_nogrp.gt.1 ) then
          CALL mp_sum(fion,3*maxsys%nax*maxsys%nsx,parai%cp_inter_grp)
       END IF
    END IF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE newd
  ! ==================================================================
  SUBROUTINE prep_bigmem_newd(nst,ig_start,nhg,na,deeq,f,fion,vpot,fnl_p,tfor)
    !driver routine for bigmem case
    !handles memory allocations and offset counters
    REAL(REAL_8),INTENT(OUT)                 :: deeq(ions1%nat,maxsys%nhxs,*)
    COMPLEX(real_8),INTENT(IN)               :: vpot(ncpw%nhg)
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: fion(:,:,:)
    INTEGER,INTENT(IN)                       :: nst(2), ig_start, nhg, na(2,ions1%nsp)
    REAL(real_8),INTENT(IN)                  :: f(*), fnl_p(il_fnl_packed(1),*)
    LOGICAL,INTENT(IN)                       :: tfor

#ifdef _USE_SCRATCHLIBRARY
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: qg1(:,:,:), vtmp(:)
    REAL(real_8), POINTER __CONTIGUOUS       :: fnlt(:,:,:),gk_trans(:,:,:),ylm(:,:)
#else
    REAL(real_8), ALLOCATABLE                :: fnlt(:,:,:),gk_trans(:,:,:),ylm(:,:)
    COMPLEX(real_8), ALLOCATABLE             :: qg1(:,:,:), vtmp(:)
#endif

    INTEGER                                  :: i, is, ia, blocksize, blocks, last_block, &
                                                isa0, nhh0, nhh, num_orb, &
                                                offset_fnl0, il_qg1(3), il_ylm(2),&
                                                il_fnlt(3), il_gk_trans(3), ierr
    CHARACTER(*), PARAMETER                  :: procedureN='prep_bigmem'

    !blocking parameters
    blocksize=cnti%blocksize_uspp*parai%ncpus
    IF(tfor)blocksize=blocksize/4
    IF(blocksize.GT.nhg)blocksize=nhg
    blocks=nhg/blocksize
    last_block=MOD(nhg,blocksize)
    !
    num_orb=nst(2)-nst(1)+1

    !Memory allocation
    il_qg1(1)=blocksize
    il_qg1(2)=0
    DO is=1,ions1%nsp
       IF(pslo_com%tvan(is))THEN
          il_qg1(2)=MAX(il_qg1(2),ions0%na(is)*(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2)
       END IF
    END DO
    il_qg1(3)=1
    il_ylm(1)=il_qg1(2)
    il_ylm(2)=1
    !if ionic forces are needed we combine both glosums:
    !ylm(:,:,1) contains rho_nm
    !ylm(:,:,2) contains local part of deeq
    !ylm(:,:,3-5) contains derivative of qg
    !ylm(:,:,2-5) enter dgemm
    !SDF allocate gk even if not requiered
    il_gk_trans(1)=blocksize
    il_gk_trans(2)=3
    il_gk_trans(3)=1
    il_fnlt(1)=1
    il_fnlt(2)=1
    il_fnlt(3)=1
    IF (tfor) THEN
       il_qg1(3)=4
       il_ylm(2)=5
       il_fnlt(1)=num_orb
       il_fnlt(2)=maxsys%nhxs
       il_fnlt(3)=parai%ncpus
       il_gk_trans(3)=blocks
       IF(last_block.GT.0) il_gk_trans(3)=blocks+1
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_fnlt,fnlt,procedureN//'_fnlt')
    CALL request_scratch(il_gk_trans,gk_trans,procedureN//'_gk_trans')
    CALL request_scratch(il_qg1,qg1,procedureN//'_qg1')
    CALL request_scratch(il_ylm,ylm,procedureN//'_ylm')
#else
    ALLOCATE(fnlt(il_fnlt(1),il_fnlt(2),il_fnlt(3)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlt)',&
         __LINE__,__FILE__)
    ALLOCATE(gk_trans(il_gk_trans(1),il_gk_trans(2), il_gk_trans(3)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate gk_trans)',&
         __LINE__,__FILE__)
    ALLOCATE(qg1(il_qg1(1),il_qg1(2),il_qg1(3)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate qg1)',&
         __LINE__,__FILE__)
    ALLOCATE(ylm(il_ylm(1),il_ylm(2)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate ylm)',&
         __LINE__,__FILE__)
#endif
    IF(tfor) CALL transpose_gk(blocks,blocksize,last_block,ig_start,gk,gk_trans)

    isa0=0
    nhh0=1
    offset_fnl0=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          nhh=(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2
          CALL evaluate_bigmem_newd(isa0,ions0%na(is),na(2,is)-na(1,is)+1,nlps_com%ngh(is),&
               nhh,offset_fnl0,na(1,is)-1,num_orb,ig_start,blocksize,last_block,blocks,tfor, &
               f(nst(1)),vpot,deeq,qg(1,nhh0),gk_trans,fnlt,qg1,ylm,&
               fnl_p(1,nst(1)),fion(:,:,is))
          nhh0=nhh0+nhh
       END IF
       isa0=isa0+ions0%na(is)
       offset_fnl0=offset_fnl0+(na(2,is)-na(1,is)+1)*nlps_com%ngh(is)
    END DO
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_ylm,ylm,procedureN//'_ylm')
    CALL free_scratch(il_qg1,qg1,procedureN//'_qg1')
    CALL free_scratch(il_gk_trans,gk_trans,procedureN//'_gk_trans')
    CALL free_scratch(il_fnlt,fnlt,procedureN//'_fnlt')
#else
    DEALLOCATE(fnlt, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fnlt)',&
         __LINE__,__FILE__)
    DEALLOCATE(gk_trans, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate gk_trans)',&
         __LINE__,__FILE__)
    DEALLOCATE(qg1, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate qg1)',&
         __LINE__,__FILE__)
    DEALLOCATE(ylm, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate ylm)',&
         __LINE__,__FILE__)
#endif

  END SUBROUTINE prep_bigmem_newd
  ! ==================================================================
  SUBROUTINE evaluate_bigmem_newd(isa0,ia_sum,ia_fnl,ngh,nhh,offset_fnl0,offset_ylm,num_orb, &
       ig_start,blocksize,last_block,blocks,tfor,f,vpot,deeq,qg_,gk_trans,fnlt,qg1,ylm, &
       fnl_p,fion)
    INTEGER,INTENT(IN)                       :: isa0, ia_sum, ia_fnl, ngh, nhh, offset_fnl0,&
                                                 offset_ylm, num_orb, ig_start, blocksize, &
                                                 last_block, blocks
    LOGICAL,INTENT(IN)                       :: tfor
    REAL(real_8),INTENT(IN)                  :: gk_trans(blocksize,3,*), f(*), &
                                                fnl_p(il_fnl_packed(1),*)
    REAL(real_8), INTENT(OUT)                :: fnlt(num_orb,ngh,*), &
                                                ylm(ia_sum,nhh,*)
    REAL(real_8), INTENT(INOUT)              :: deeq(ions1%nat,maxsys%nhxs,*)
    REAL(real_8), INTENT(INOUT) __CONTIGUOUS :: fion(:,:)
    COMPLEX(real_8), INTENT(OUT)             :: qg1(blocksize,nhh,*)
    COMPLEX(real_8), INTENT(IN)              :: qg_(ncpw%nhg,*), vpot(*)


    INTEGER                                  :: istart, iblock, qgstart, ijv, methread, ia, &
                                                isa, iv, jv, k, len_qg1, start_ylm
    REAL(real_8)                             :: fac, omtpiba, fiont(3), otr

    len_qg1=1
    start_ylm=1
    IF(tfor)THEN
       len_qg1=4
       start_ylm=2
    END IF
    omtpiba=parm%omega*parm%tpiba
    istart=ig_start-1
    qgstart=ig_start-1
    fac=0.0_real_8
    DO iblock=1,blocks
       CALL evaluate_block_newd(blocksize,blocksize,istart,qgstart,isa0,ia_sum,start_ylm, &
            nhh,fac,tfor,gk_trans(1,1,iblock),qg1,ylm,vpot(istart+1),qg_(1,1))

       IF(iblock.EQ.1)THEN
          fac=1.0_real_8
          !remove double counting of geq0
          IF(geq0.AND.istart.EQ.0) THEN
             CALL dger(ia_sum,len_qg1*nhh,-1.0_real_8,&
                  eigrb(1,isa0+1),2*ncpw%nhg, &
                  qg1,2*blocksize,&
                  ylm(1,1,start_ylm),ia_sum)

          END IF
       END IF
       istart=istart+blocksize
       qgstart=qgstart+blocksize
    END DO
    IF(last_block.GT.0)THEN
       iblock=blocks+1
       CALL evaluate_block_newd(blocksize,last_block,istart,qgstart,isa0,ia_sum,start_ylm, &
            nhh,fac,tfor,gk_trans(1,1,iblock),qg1,ylm,vpot(istart+1),qg_(1,1))
    END IF

    IF(tfor)THEN
       ijv=0
       methread=1
       !$omp parallel private(methread)
       !$ methread=omp_get_thread_num()+1
       CALL calc_rho(ngh,ia_fnl,offset_fnl0,num_orb,offset_ylm,ia_sum,&
            ylm,fnlt(:,:,methread),f,fnl_p)
       !$omp end parallel
    ENDIF
    CALL mp_sum(ylm,nhh*ia_sum*start_ylm,parai%cp_grp)
    !$omp parallel private(ijv,iv,jv,ia,isa)
    ijv=0
    DO iv=1,ngh
       DO jv=iv,ngh
          ijv=ijv+1
          !$omp do
          DO ia=1,ia_sum
             isa=isa0+ia
             deeq(isa,iv,jv)=ylm(ia,ijv,start_ylm)*parm%omega
             deeq(isa,jv,iv)=ylm(ia,ijv,start_ylm)*parm%omega
          ENDDO
          !$omp end do nowait
       ENDDO
    ENDDO
    !$omp end parallel
    IF (tfor) THEN
       !$omp parallel do private(ia,fiont,ijv,otr,k)
       DO ia=1,ia_sum
          fiont=0._real_8
          ijv=0
          DO ijv=1,nhh
             OTR=YLM(ia,ijv,1)*omtpiba
             DO k=1,3
                fiont(k)=fiont(k)+OTR*YLM(ia,ijv,2+k)
             ENDDO
          END DO
          DO k=1,3
             fion(k,ia)=fion(k,ia)+fiont(k)
          END DO
       END DO
    ENDIF
  END SUBROUTINE evaluate_bigmem_newd

  SUBROUTINE transpose_gk(blocks,blocksize,last_block,ig_start,gk_,gktrans)
    ! ==---------------------------------------------------------_-----==
    INTEGER,INTENT(IN)                         :: blocks, blocksize, last_block, ig_start
    REAL(real_8),INTENT(IN)                    :: gk_(3,ncpw%nhg)
    REAL(real_8),INTENT(OUT)                   :: gktrans(blocksize,3,*)

    INTEGER                                    :: istart,iblock,ig,k
    !$omp parallel private(istart,iblock,ig,k)
    istart=ig_start-1
    DO iblock=1,blocks
       !$omp do __COLLAPSE2
       DO ig=1,blocksize
          DO k=1,3
             gktrans(ig,k,iblock)=gk(k,istart+ig)
          END DO
       END DO
       !$omp end do nowait
       istart=istart+blocksize
    END DO
    IF(last_block.GT.0)THEN
       iblock=blocks+1
       !$omp do __COLLAPSE2
       DO ig=1,last_block
          DO k=1,3
             gktrans(ig,k,iblock)=gk(k,istart+ig)
          END DO
       END DO
       !$omp end do nowait
    END IF
    !$omp end parallel
  END SUBROUTINE transpose_gk
  ! ==================================================================
  SUBROUTINE calc_rho(ngh,ia_fnl,offset_fnl0,norb,offset_ylm,ld_ylm,ylm,temp,f,fnl_p)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                         :: ngh,ia_fnl,offset_fnl0,norb,offset_ylm,&
                                                  ld_ylm
    REAL(real_8),INTENT(IN)                    :: fnl_p(il_fnl_packed(1),norb),f(norb)
    REAL(real_8),INTENT(OUT)                   :: ylm(ld_ylm,*),temp(norb,*)
    INTEGER                                    :: fnl_offset,ia,iv,i,ijv,jv,fnl_ia
    REAL(real_8)                               :: ftmp,tmp

    !$omp do
    DO ijv=1,ngh*(ngh+1)/2
       DO ia=1,ld_ylm
          ylm(ia,ijv)=0.0_real_8
       END DO
    END DO
    !$omp do
    DO ia=1,ia_fnl
       fnl_offset=offset_fnl0
       DO iv=1,ngh
          fnl_ia=ia+fnl_offset
          DO i=1,norb
             temp(i,iv)=fnl_p(fnl_ia,i)
          END DO
          fnl_offset=fnl_offset+ia_fnl
       END DO
       ijv=0
       DO iv=1,ngh
          DO jv=iv,ngh
             ijv=ijv+1
             tmp=0._real_8
             DO i=1,norb
                tmp=tmp+f(i)*temp(i,iv)*temp(i,jv)
             END DO
             ftmp=1.0_real_8
             IF(iv.NE.jv) ftmp=2._real_8
             YLM(ia+offset_ylm,ijv)=tmp*ftmp
          END DO
       END DO
    END DO
  END SUBROUTINE calc_rho
  ! ==================================================================
  SUBROUTINE evaluate_block_newd(ld_qg1,blocksize,block_start,qgstart,isa0,ia_sum,start_ylm, &
       nhh,fac,tfor,gk_trans,qg1,ylm,vpot,qg_)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                         :: ld_qg1, blocksize, block_start, isa0, &
                                                  ia_sum, start_ylm, nhh, qgstart
    LOGICAL,INTENT(IN)                         :: tfor
    REAL(real_8),INTENT(IN)                    :: fac, gk_trans(ld_qg1,3)
    COMPLEX(real_8),INTENT(IN)                 :: vpot(blocksize), qg_(ncpw%nhg,*)
    COMPLEX(real_8),INTENT(OUT)                :: qg1(blocksize,nhh,*)
    REAL(real_8)                               :: ylm(ia_sum,nhh,*)
    INTEGER                                    :: ijv,ig,k,ig2,len_qg1

    !$omp parallel private(ijv,ig,k,ig2)
    DO ijv=1,nhh
       IF(tfor)THEN
          !$omp do
          DO ig=1,blocksize
             ig2=qgstart+ig
             qg1(ig,ijv,1)=CONJG(qg_(ig2,ijv))*vpot(ig)
             DO k=1,3
                qg1(ig,ijv,1+k)=&
                     CMPLX(gk_trans(ig,k)*AIMAG(qg1(ig,ijv,1)),&
                     gk_trans(ig,k)*(-REAL(qg1(ig,ijv,1),kind=real_8)))
             END DO
          END DO
          !$omp end do nowait
       ELSE
          !$omp do
          DO ig=1,blocksize
             ig2=qgstart+ig
             qg1(ig,ijv,1)=CONJG(qg_(ig2,ijv))*vpot(ig)
          ENDDO
          !$omp end do nowait
       END IF
    ENDDO
    !$omp end parallel
    len_qg1=1
    IF(tfor)len_qg1=4
    CALL dgemm('T','N',ia_sum,nhh*len_qg1,2*blocksize,2.0_real_8,&
         eigrb(block_start+1,isa0+1),2*ncpw%nhg,&
         qg1,2*blocksize,&
         fac,ylm(1,1,start_ylm),ia_sum)
    RETURN
  END SUBROUTINE evaluate_block_newd
  ! ==================================================================
  SUBROUTINE prep_smallmem_newd(nst,ig_start,nhg,na,deeq,f,fion,vpot,fnl_p,tfor)
    REAL(REAL_8),INTENT(OUT)                 :: deeq(ions1%nat,maxsys%nhxs,*)
    COMPLEX(real_8),INTENT(IN)               :: vpot(ncpw%nhg)
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: fion(:,:,:)
    INTEGER,INTENT(IN)                       :: nst(2), na(2,ions1%nsp), ig_start, nhg
    REAL(real_8),INTENT(IN)                  :: f(*), fnl_p(il_fnl_packed(1),*)
    LOGICAL,INTENT(IN)                       :: tfor

    COMPLEX(real_8), ALLOCATABLE             :: vtmp(:)
    REAL(real_8), ALLOCATABLE                :: ylm(:), fnlt(:,:,:)
    INTEGER                                  :: isa0, offset_fnl0, is, ia_sum, ijv, iv, jv, &
                                                offset_iv, offset_jv, ig, ia, isa, i, k, &
                                                nhh, methread, offset_ylm, il_vtmp(1), &
                                                il_ylm(1), il_fnlt(3), ierr
    CHARACTER(*), PARAMETER                  :: procedureN='smallmem_newd'
    REAL(real_8)                             :: ftmp, otr
    REAL(real_8), EXTERNAL                   :: dotp
    LOGICAL                                  :: geq0_back
    ! ==--------------------------------------------------------------==

    il_ylm(1)=0
    IF(tfor)THEN
       DO is=1,ions1%nsp
          IF(pslo_com%tvan(is))THEN
             il_ylm(1)=MAX(il_ylm(1),ions0%na(is)*(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2)
          END IF
       END DO
    END IF
    il_vtmp(1)=nhg
    IF(tfor)THEN
       il_fnlt(1)=nst(2)-nst(1)+1
       il_fnlt(2)=maxsys%nhxs
       il_fnlt(3)=parai%ncpus
       ALLOCATE(fnlt(il_fnlt(1),il_fnlt(2),il_fnlt(3)), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlt)',&
            __LINE__,__FILE__)
    END IF

    ALLOCATE(qg(ncpw%nhg,1), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate qg)',&
         __LINE__,__FILE__)
    ALLOCATE(vtmp(il_vtmp(1)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate vtmp)',&
         __LINE__,__FILE__)
    ALLOCATE(ylm(il_ylm(1)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate ylm)',&
         __LINE__,__FILE__)

    !dotp function is not cp_grp aware... temporary redefine geq0
    geq0_back=geq0
    geq0=(ig_start.EQ.1.AND.geq0)
    !
    isa0=0
    offset_fnl0=0
    DO is=1,ions1%nsp
       ia_sum=na(2,is)-na(1,is)+1
       IF (pslo_com%tvan(is)) THEN
          IF(tfor)THEN
             nhh=(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2
             methread=1
             !$omp parallel private(methread)
             !$ methread=omp_get_thread_num()+1
             CALL calc_rho(nlps_com%ngh(is),ia_sum,offset_fnl0,nst(2)-nst(1)+1,&
                  na(1,is)-1,ions0%na(is),ylm,fnlt(:,:,methread),f(nst(1)),&
                  fnl_p(1,nst(1)))
             !$omp end parallel
             CALL mp_sum(ylm,nhh*ions0%na(is),parai%cp_grp)
          END IF

          ijv=0
          DO iv=1,nlps_com%ngh(is)
             DO jv=iv,nlps_com%ngh(is)
                ijv=ijv+1
                CALL qvan2(iv,jv,is,qg(:,1))
                !$omp parallel do private(IG) shared(QG) schedule(static)
                DO ig=ig_start,ig_start+nhg-1
                   qg(ig,1)=CONJG(qg(ig,1))*vpot(ig)
                ENDDO
                DO ia=1,ions0%na(is)
                   isa=isa0+ia
                   !$omp parallel do private(IG) shared(VTMP) schedule(static)
                   DO ig=ig_start,ig_start+nhg-1
                      vtmp(ig-ig_start+1)=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                           ei3(isa,inyh(3,ig))
                   ENDDO
                   deeq(isa,iv,jv)=parm%omega*dotp(nhg,qg(ig_start,1),vtmp)
                   deeq(isa,jv,iv)=deeq(isa,iv,jv)
                   IF (tfor) THEN
                      offset_ylm=ions0%na(is)*(ijv-1)+ia
                      otr=parm%omega*parm%tpiba*ylm(offset_ylm)
                      DO k=1,3
                         CALL cftemp(nhg,qg(ig_start,1),vtmp,gk(1,ig_start),k,ftmp)
!                         IF (iv.NE.jv) ftmp=2.0_real_8*ftmp
                         fion(k,ia,is)=fion(k,ia,is)+ftmp*otr
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       offset_fnl0=offset_fnl0+(na(2,is)-na(1,is)+1)*nlps_com%ngh(is)
       isa0=isa0+ions0%na(is)
    ENDDO
    geq0=geq0_back
    DEALLOCATE(qg, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate qg)',&
         __LINE__,__FILE__)
    DEALLOCATE(vtmp, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate vtmp)',&
         __LINE__,__FILE__)
    DEALLOCATE(ylm, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate ylm)',&
         __LINE__,__FILE__)
    IF(tfor)THEN
       DEALLOCATE(fnlt, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlt)',&
            __LINE__,__FILE__)
    END IF
    CALL mp_sum(deeq,ions1%nat*maxsys%nhxs*maxsys%nhxs,parai%cp_grp)
  END SUBROUTINE prep_smallmem_newd
END MODULE newd_utils

! ==================================================================
SUBROUTINE cftemp(nhg,qg,vtmp,gk,k,ftmp)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: nhg
  REAL(real_8)                               :: qg(2,nhg), vtmp(2,nhg), &
                                                gk(3,nhg)
  INTEGER                                    :: k
  REAL(real_8)                               :: ftmp

  INTEGER                                    :: ig

! Variables
! ==--------------------------------------------------------------==
! G=0 TERM IS ZERO !

  ftmp=0.0_real_8
  !$omp parallel do private(IG) reduction(+:FTMP) schedule(static)
  DO ig=1,nhg
     ftmp=ftmp+gk(k,ig)*(vtmp(1,ig)*qg(2,ig)-vtmp(2,ig)*qg(1,ig))
  ENDDO
  ftmp=2.0_real_8*ftmp
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE cftemp
