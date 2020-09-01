#include "cpmd_global.h"

MODULE rhov_utils
  USE cppt,                            ONLY: indz,&
                                             inyh,&
                                             nzh
  USE cp_grp_utils,                    ONLY: cp_grp_redist_array,&
                                             cp_grp_get_sizes,&
                                             cp_grp_split_atoms
  USE cvan,                            ONLY: qg
  USE distribution_utils,              ONLY: dist_entity
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE fft_maxfft,                      ONLY: maxfftn
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nlps_com
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
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw, &
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

  PUBLIC :: rhov

CONTAINS

  ! ==================================================================
  SUBROUTINE rhov(i_start,i_end,rsumv,psi)
    ! ==--------------------------------------------------------------==
    ! ==         THIS ROUTINE CALCULATES THE VANDERBILT DENSITY       ==
    ! ==                                                              ==
    ! == N_V(G) = SUM_I,IJ RHO_I,IJ Q_I,JI(G) E^-IG.R_I               ==
    ! == RHO_I,IJ = SUM_N < BETA_I,I | PSI_N >< PSI_N | BETA_I,J >    ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: i_start, i_end
    REAL(real_8),INTENT(OUT)                 :: rsumv
    COMPLEX(real_8),INTENT(INOUT), OPTIONAL &
         __CONTIGUOUS                        :: psi(:)

    INTEGER                                  :: ig, isub, isub1, ig_start, nhg_loc, &
                                                ierr, isub2
    INTEGER, ALLOCATABLE                     :: na_grp(:,:,:), na(:,:), nst(:,:)
    INTEGER(int_8)                           :: il_deltar(1)
#ifdef _USE_SCRATCHLIBRARY
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: deltar(:)
#else
    COMPLEX(real_8), ALLOCATABLE             :: deltar(:)
#endif
    CHARACTER(*), PARAMETER                  :: procedureN='rhov'
    IF(cntl%fft_tune_batchsize)THEN
       CALL tiset(procedureN//'tuning',isub2)
    ELSE
       CALL tiset(procedureN,isub)
    END IF

    IF (cntl%tfdist) CALL stopgm(procedureN,'TFDIST NOT IMPLEMENTED',&
         __LINE__,__FILE__)
    IF (imagp.EQ.2)CALL stopgm(procedureN,'K-POINT NOT IMPLEMENTED',&
         __LINE__,__FILE__)
    CALL setfftn(0)
    ALLOCATE(na_grp(2,ions1%nsp,0:parai%cp_nogrp-1), na(2,ions1%nsp), nst(2,0:parai%nproc-1),&
         stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'allocation problem)',&
         __LINE__,__FILE__)

    CALL cp_grp_get_sizes(first_nhg=ig_start,nhg_l=nhg_loc)
    CALL dist_entity(i_end-i_start+1,parai%nproc,nst)
    !shift everything to starting state
    nst(1,parai%me)=nst(1,parai%me)+i_start-1
    nst(2,parai%me)=nst(2,parai%me)+i_start-1
    CALL cp_grp_split_atoms(na_grp)
    na(:,:)=na_grp(:,:,parai%cp_inter_me)


    il_deltar(1)=ncpw%nhg
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_deltar,deltar,procedureN//'_deltar',ierr)
#else
    ALLOCATE(deltar(il_deltar(1)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
   ! ==--------------------------------------------------------------==
    IF (cntl%bigmem) THEN
       CALL prep_bigmem_rhov(nst(1,parai%me),ig_start,nhg_loc,na,deltar,fnl_packed)
    ELSE
       CALL prep_smallmem_rhov(nst(1,parai%me),ig_start,nhg_loc,na,deltar,fnl_packed)
    ENDIF
    IF (parai%cp_nogrp.GT.1) THEN
       CALL TISET(procedureN//'_grpsb',isub1)
       CALL cp_grp_redist_array(deltar,ncpw%nhg)
       CALL TIHALT(procedureN//'_grpsb',isub1)
    END IF
    IF (geq0) THEN
       rsumv=REAL(deltar(1))
    ELSE
       rsumv=0.0_real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (PRESENT(psi)) THEN
       CALL zeroing(psi)!,maxfft)
       !$omp parallel do private(IG) shared(PSI)
       DO ig=1,ncpw%nhg
          psi(nzh(ig))=deltar(ig)
          psi(indz(ig))=CONJG(deltar(ig))
       ENDDO
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_deltar,deltar,procedureN//'_deltar',ierr)
#else
    DEALLOCATE(deltar,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    IF (PRESENT(psi)) THEN
       CALL invfftn(psi, .FALSE.,parai%allgrp)
    END IF
    DEALLOCATE(na_grp, na, nst, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'deallocation problem)',&
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
    IF(cntl%fft_tune_batchsize)THEN
       CALL tihalt(procedureN//'tuning',isub2)
    ELSE
       CALL tihalt(procedureN,isub)
    END IF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rhov
  ! ==================================================================
  SUBROUTINE prep_bigmem_rhov(nst,ig_start,nhg,na,deltar,fnl_p)
    !driver routine for bigmem case
    !handles memory allocations and offset counters
    COMPLEX(real_8),INTENT(OUT)              :: deltar(ncpw%nhg)
    INTEGER,INTENT(IN)                       :: nst(2), ig_start, nhg, na(2,ions1%nsp)
    REAL(real_8),INTENT(IN)                  :: fnl_p(il_fnl_packed(1),*)

    INTEGER                                  :: is, isa0, nhh0, offset_fnl0, nhh, num_orb, &
                                                blocksize, blocks, last_block, ig, ierr
    INTEGER(int_8)                           :: il_dia(1), il_ctmp(2), il_fnlt(3)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: dia(:), fnlt(:,:,:)
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: ctmp(:,:)
#else
    REAL(real_8), ALLOCATABLE                :: dia(:), fnlt(:,:,:)
    COMPLEX(real_8), ALLOCATABLE             :: ctmp(:,:)
#endif
    CHARACTER(*), PARAMETER                  :: procedureN='prep_bigmem_rhov'


    blocksize=cnti%blocksize_uspp*parai%ncpus
    IF(blocksize.GT.nhg)blocksize=nhg
    blocks=nhg/blocksize
    last_block=MOD(nhg,blocksize)

    num_orb=nst(2)-nst(1)+1

    il_dia(1)=0
    DO is=1,ions1%nsp
       IF(pslo_com%tvan(is))THEN
          il_dia(1)=MAX(il_dia(1),ions0%na(is)*(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2)
       END IF
    END DO

    il_ctmp(1)=blocksize
    il_ctmp(2)=(maxsys%nhxs*(maxsys%nhxs+1))/2
    il_fnlt(1)=num_orb
    il_fnlt(2)=maxsys%nhxs
    il_fnlt(3)=parai%ncpus
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_dia,dia,procedureN//'_dia',ierr)
#else
    ALLOCATE(dia(il_dia(1)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_fnlt,fnlt,procedureN//'_fnlt',ierr)
#else
    ALLOCATE(fnlt(il_fnlt(1),il_fnlt(2),il_fnlt(3)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_ctmp,ctmp,procedureN//'_ctmp',ierr)
#else
    ALLOCATE(ctmp(il_ctmp(1),il_ctmp(2)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    !$omp parallel do private(ig)
    DO ig=ig_start,ig_start+nhg-1
       deltar(ig)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
    END DO

    isa0=0
    nhh0=1
    offset_fnl0=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          nhh=(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2
          CALL evaluate_bigmem_rhov(isa0,ions0%na(is),na(2,is)-na(1,is)+1,nlps_com%ngh(is),&
               nhh,offset_fnl0,na(1,is)-1,num_orb,ig_start,blocksize,last_block,blocks,&
               crge%f(nst(1),1),dia,deltar,qg(1,nhh0),fnlt,fnl_p(1,nst(1)),ctmp)
          nhh0=nhh0+nhh
       END IF
       isa0=isa0+ions0%na(is)
       offset_fnl0=offset_fnl0+(na(2,is)-na(1,is)+1)*nlps_com%ngh(is)
    END DO

#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_ctmp,ctmp,procedureN//'_ctmp',ierr)
#else
    DEALLOCATE(ctmp,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_fnlt,fnlt,procedureN//'_fnlt',ierr)
#else
    DEALLOCATE(fnlt,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_dia,dia,procedureN//'_dia',ierr)
#else
    DEALLOCATE(dia,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

  END SUBROUTINE prep_bigmem_rhov
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE evaluate_bigmem_rhov(isa0,ia_sum,ia_fnl,ngh,nhh,offset_fnl0,offset_ylm,num_orb,&
       ig_start,blocksize,last_block,blocks,f,dia,deltar,qg_,fnlt,fnl_p,ctmp)
    INTEGER,INTENT(IN)                       :: isa0, ia_sum, ia_fnl, ngh, nhh, offset_fnl0,&
                                                offset_ylm, num_orb, ig_start, blocksize, &
                                                last_block, blocks
    REAL(real_8), INTENT(OUT)                :: dia(ia_sum,nhh), fnlt(num_orb,ngh,*)
    REAL(real_8), INTENT(IN)                 :: fnl_p(il_fnl_packed(1),num_orb), f(num_orb)
    COMPLEX(real_8),INTENT(IN)               :: qg_(ncpw%nhg,*)
    COMPLEX(real_8),INTENT(OUT)              :: ctmp(blocksize,nhh)
    COMPLEX(real_8),INTENT(INOUT)            :: deltar(*)

    INTEGER                                  :: methread, istart, qgstart, iblock
    methread=1
    !$omp parallel private(methread)
    !$ methread=omp_get_thread_num()+1
    CALL calc_rho(ngh,ia_fnl,offset_fnl0,num_orb,offset_ylm,ia_sum,dia,fnlt(:,:,methread),&
         f,fnl_p)
    !$omp end parallel
    CALL mp_sum(dia,nhh*ia_sum,parai%cp_grp)
    istart=ig_start
    qgstart=ig_start-1
    DO iblock=1,blocks
       CALL evaluate_block_rhov(blocksize,istart,qgstart,isa0,ia_sum,nhh,dia,ctmp,&
            deltar(istart),qg_)
       qgstart=qgstart+blocksize
       istart=istart+blocksize
    END DO
    IF(last_block.gt.0)THEN
       CALL evaluate_block_rhov(last_block,istart,qgstart,isa0,ia_sum,nhh,dia,ctmp,&
            deltar(istart),qg_)
    END IF


  END SUBROUTINE evaluate_bigmem_rhov
  ! ==================================================================
  SUBROUTINE evaluate_block_rhov(blocksize,block_start,qgstart,isa0,ia_sum,nhh,dia,ctmp,&
       deltar,qg_)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                         :: blocksize, block_start, qgstart, isa0, &
                                                  ia_sum, nhh
    REAL(real_8),INTENT(IN)                    :: dia(ia_sum,nhh)
    COMPLEX(real_8),INTENT(IN)                 :: qg_(ncpw%nhg,*)
    COMPLEX(real_8),INTENT(OUT)                :: ctmp(blocksize,nhh)
    COMPLEX(real_8),INTENT(INOUT)              :: deltar(blocksize)
    INTEGER                                    :: ijv, ig, ig2

    CALL dgemm('n','n',2*blocksize,nhh,ia_sum &
         ,1.0d0,EIGRB(block_start,isa0+1),2*ncpw%nhg &
         ,dia,ia_sum,0.0_real_8,CTMP,2*blocksize)
    !$omp parallel private(ijv,ig,ig2)
    DO ijv=1,nhh
       !$omp do
       DO ig=1,blocksize
          ig2=qgstart+ig
          deltar(ig)=deltar(ig)+qg_(ig2,ijv)*ctmp(ig,ijv)
       ENDDO
       !$omp end do nowait
    END DO
    !$omp end parallel
    RETURN
  END SUBROUTINE evaluate_block_rhov
  ! ==================================================================
  SUBROUTINE prep_smallmem_rhov(nst,ig_start,nhg,na,deltar,fnl_p)
    COMPLEX(real_8),INTENT(OUT)              :: deltar(ncpw%nhg)
    INTEGER,INTENT(IN)                       :: nst(2), ig_start, nhg, na(2,ions1%nsp)
    REAL(real_8),INTENT(IN)                  :: fnl_p(il_fnl_packed(1),*)

    INTEGER                                  :: isa0, is, iv, jv, ia, isa, i, ig, nhh, ijv, &
                                                methread, offset_fnl0, offset_dia, &
                                                il_ctmp(1), il_dia(1), il_fnlt(3), ierr
    REAL(real_8)                             :: sum
    REAL(real_8), ALLOCATABLE                :: dia(:), fnlt(:,:,:)
    COMPLEX(real_8), ALLOCATABLE             :: ctmp(:)
    CHARACTER(*), PARAMETER                  :: procedureN='prep_smallmem_rhov'


    il_dia(1)=0
    DO is=1,ions1%nsp
       IF(pslo_com%tvan(is))THEN
          il_dia(1)=MAX(il_dia(1),ions0%na(is)*(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2)
       END IF
    END DO
    il_fnlt(1)=nst(2)-nst(1)+1
    il_fnlt(2)=maxsys%nhxs
    il_fnlt(3)=parai%ncpus
    il_ctmp(1)=ncpw%nhg

    ALLOCATE(dia(il_dia(1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(fnlt(il_fnlt(1),il_fnlt(2),il_fnlt(3)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ctmp(il_ctmp(1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(qg(ncpw%nhg,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    !$omp parallel do private(ig)
    DO ig=ig_start,ig_start+nhg-1
       deltar(ig)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
    END DO

    isa0=0
    offset_fnl0=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          nhh=(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2
          methread=1
          !$omp parallel private(methread)
          !$ methread=omp_get_thread_num()+1
          CALL calc_rho(nlps_com%ngh(is),na(2,is)-na(1,is)+1,offset_fnl0,nst(2)-nst(1)+1, &
               na(1,is)-1,ions0%na(is),dia,fnlt(:,:,methread),crge%f(nst(1),1),&
               fnl_p(1,nst(1)))
          !$omp end parallel
          CALL mp_sum(dia,nhh*ions0%na(is),parai%cp_grp)
          ijv=0
          DO iv=1,nlps_com%ngh(is)
             DO jv=iv,nlps_com%ngh(is)
                ijv=ijv+1
                CALL qvan2(iv,jv,is,qg(:,1))
                DO ia=1,ions0%na(is)
                   isa=isa0+ia
                   !$omp parallel do private(IG)
                   DO ig=ig_start,ig_start+nhg-1
                      ctmp(ig)=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                           ei3(isa,inyh(3,ig))
                   ENDDO
                   offset_dia=ions0%na(is)*(ijv-1)+ia
                   !$omp parallel do private(IG)
                   DO ig=ig_start,ig_start+nhg-1
                      deltar(ig)=deltar(ig)+qg(ig,1)*dia(offset_dia)*ctmp(ig)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       offset_fnl0=offset_fnl0+(na(2,is)-na(1,is)+1)*nlps_com%ngh(is)
       isa0=isa0+ions0%na(is)
    ENDDO
    DEALLOCATE(qg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ctmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(dia,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(fnlt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

  END SUBROUTINE prep_smallmem_rhov
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

END MODULE rhov_utils
