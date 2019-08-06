#include "cpmd_global.h"

MODULE rnlsm2_utils
  USE cppt,                            ONLY: gk,&
                                             twnl
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms,&
                                             cp_grp_redist_dfnl_fnl
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE fnl_utils,                       ONLY: unpack_dfnl,&
                                             sort_dfnl
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8,&
                                             int_8,&
                                             int_4
  USE kpnt,                            ONLY: eigkr,&
                                             rk
  USE kpts,                            ONLY: tkpts
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_sum,mp_bcast
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlm
  !$ USE omp_lib,                      ONLY: omp_set_nested,&
  !$                                         omp_get_thread_num,&
  !$                                         omp_set_num_threads
  USE parac,                           ONLY: parai,&
                                             paral
  USE reshaper,                        ONLY: reshape_inplace
  USE rnlsm_helper,                    ONLY: proj_beta,&
                                             set_buffers,&
                                             tune_rnlsm
  USE sfac,                            ONLY: dfnl,&
                                             dfnla,&
                                             dfnl_packed,&
                                             eigr,&
                                             il_dfnl_packed
  USE system,                          ONLY: cntl,&
                                             cnti,&
                                             cntr,&
                                             ncpw,&
                                             nkpt,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlsm2

CONTAINS

  ! ==================================================================
  SUBROUTINE RNLSM2(c0,nstate,ikpt,ikind,dfnl_unpack)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ARRAY DFNL WHICH IS USED IN THE SUBROUTINE RNLFOR       ==
    ! ==  K-POINT VERSION IS IMPLEMENTED                              ==
    ! ==          NOT IMPLEMENTED FOR TSHEL(IS)=TRUE                  ==
    ! == Rewritten: Tobias Kloeffel, FAU Erlangen-Nuernberg April 2019==
    ! == Introduce overlapping communication computation algorithm    ==
    ! == enable via USE_OVERLAPPING_COMM_COMP in the &CPMD section    ==
    ! == Set block counts via INPUT: RNLSM2_BLOCKCOUNT                ==
    ! == First blocksize: RNLSM2_BLOCKSIZE1 (non overlapping part)    ==
    ! == Other blocksizes: RNLSM2_BLOCKSIZE2 (overlapping part)       ==
    ! == Or use autotuning: RNLSM2_AUTOTUNE (needs number of          ==
    ! == iterations for timing measurements                           ==
    ! == special case for uspp: keep only a packed version of the dfnl==
    ! == array, no glosum performed, optimzed rnlfor/rnlfl routines   ==
    ! == are way faster than a global summation                       ==
    ! == cp_grp distribution along beta projectors, only active for   ==
    ! == uspp part(disabled unpacking of dfnl array)                  ==
    ! == full performance only with saved arrays or scratch_library   ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate, ikpt, ikind
    COMPLEX(real_8),INTENT(IN)               :: c0(nkpt%ngwk,nstate)
    LOGICAL,OPTIONAL,INTENT(IN)              :: dfnl_unpack

    LOGICAL                                  :: unpack
    CHARACTER(*), PARAMETER                  :: procedureN = 'rnlsm2'
    REAL(real_8)                             :: temp
    INTEGER, PARAMETER                       :: maxbuff=15
    INTEGER                                  :: i, ierr, ig, k, isub, isub2, &
                                                methread, nthreads, nested_threads, &
                                                tot_work, start_dai, ld_dai, end_dai, &
                                                ia_sum,start_ia,end_ia,&
                                                buffcount, buff, igeq0, ld_buffer(maxbuff),&
                                                start_buffer(maxbuff)
    INTEGER(int_8)                           :: il_gktemp(2), il_eiscr(2), il_t(1), il_dai(1)
    INTEGER,ALLOCATABLE                      :: na_buff(:,:,:), na_grp(:,:,:)
    REAL(real_8),POINTER __CONTIGUOUS        :: dai(:)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8),POINTER __CONTIGUOUS        :: t(:),gktemp(:,:)
    COMPLEX(real_8),POINTER __CONTIGUOUS     :: eiscr(:,:)
#else
    REAL(real_8),ALLOCATABLE                 :: t(:),gktemp(:,:)
    COMPLEX(real_8),ALLOCATABLE              :: eiscr(:,:)
#endif
    REAL(real_8), SAVE                       :: timings(2)=0.0_real_8
    INTEGER, SAVE                            :: autotune_it=0
    !$ LOGICAL                               :: locks(maxbuff)

    ! ==--------------------------------------------------------------==
    ! IF no non-local components -> return.
    IF(NLM.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL TISET(procedureN,ISUB)
    ! ==--------------------------------------------------------------==

    ! If we store packed dfnl, we are able to skip glosum
    IF(PRESENT(dfnl_unpack))THEN
       unpack=dfnl_unpack
    ELSE
       unpack=.TRUE.
    END IF
    ! split atoms between cp_groups
    ALLOCATE(na_grp(2,ions1%nsp,0:parai%cp_nogrp-1), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate na_grp',&
         __LINE__,__FILE__)

    CALL cp_grp_split_atoms(na_grp)

    ! set autotuning stuff
    CALL set_buffers(autotune_it,maxbuff,timings,cnti%rnlsm2_bc,cntr%rnlsm2_b1,&
         cntr%rnlsm2_b2,na_grp(:,:,parai%cp_inter_me),na_buff,ld_buffer,start_buffer,&
         tot_work,nstate,3*imagp,.NOT.unpack)
    buffcount=cnti%rnlsm2_bc
    ! allocate memory
    il_dai(1)=tot_work*nstate
    il_eiscr(1)=nkpt%ngwk
    il_eiscr(2)=MAXVAL(ld_buffer)/imagp
    il_gktemp(1)=nkpt%ngwk
    il_gktemp(2)=3
    il_t(1)=nkpt%ngwk
    IF(buffcount.GT.1)THEN
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_dai,dai,procedureN//'_dai')
#else
       ALLOCATE(dai(il_dai(1)), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate dai',&
            __LINE__,__FILE__)
#endif
    ELSE
       CALL reshape_inplace(dfnl_packed, (/INT(il_dai(1),kind=int_4)/), dai)
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_eiscr,eiscr,procedureN//'_eiscr')
    CALL request_scratch(il_gktemp,gktemp,procedureN//'_gktemp')
    CALL request_scratch(il_t,t,procedureN//'_t')
#else
    ALLOCATE(eiscr(il_eiscr(1),il_eiscr(2)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate eiscr',&
         __LINE__,__FILE__)
    ALLOCATE(gktemp(il_gktemp(1),il_gktemp(2)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate gktemp',&
         __LINE__,__FILE__)
    ALLOCATE(t(il_t(1)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate t',&
         __LINE__,__FILE__)
#endif
    !$omp parallel private(ig,k)
    IF (tkpts%tkpnt) THEN
       !$OMP do __COLLAPSE2
       DO ig=1,ncpw%ngw
          DO k=1,3
             gktemp(ig,k)=gk(k,ig)+RK(K,IKIND)
          ENDDO
       ENDDO
       !$omp end do nowait
       !$OMP do __COLLAPSE2
       DO ig=ncpw%ngw+1,nkpt%ngwk
          DO k=1,3
             gktemp(ig,k)=-gk(k,ig-ncpw%ngw)+RK(K,IKIND)
          ENDDO
       ENDDO
       !$omp end do nowait
    ELSE
       !$OMP do __COLLAPSE2
       DO ig=1,nkpt%ngwk
          DO k=1,3
             gktemp(ig,k)=gk(k,ig)
          ENDDO
       ENDDO
       !$omp end do nowait
    ENDIF
    !$omp end parallel
    IF(tkpts%tkpnt)THEN
       igeq0=ncpw%ngw+1
    ELSE
       igeq0=1
    END IF
    buff=1
    start_dai=start_buffer(buff)
    ld_dai=ld_buffer(buff)
    end_dai=start_dai-1+ld_dai*nstate
    IF (autotune_it.GT.0.AND.autotune_it.LE.cnti%rnlsm_autotune_maxit) temp=m_walltime()
    IF(ld_dai.GT.0)THEN
       CALL proj_beta(na_buff(:,:,buff),igeq0,nstate,c0,nkpt%ngwk,eigkr(:,:,ikind),&
            twnl(:,:,:,ikind),eiscr,t,nkpt%ngwk,1,dai(start_dai:end_dai),ld_dai/imagp,&
            .TRUE.,tkpts%tkpnt,geq0,gktemp)
    END IF
    IF(autotune_it.GT.0.AND.autotune_it.LE.cnti%rnlsm_autotune_maxit)&
         timings(1)=timings(1)+m_walltime()-temp

    !old algorithm active if dfnl needs to be unpacked
    IF(unpack)THEN
       IF(autotune_it.GT.0.AND.autotune_it.LE.cnti%rnlsm_autotune_maxit)&
            timings(1)=timings(1)+m_walltime()-temp

       !now we split up the threads, thread=0 is used to communicate,
       !thread=1 is used to do the dgemms and copy the data to its destionation.

       IF(cntl%overlapp_comm_comp)THEN
          nthreads=MIN(2,parai%ncpus)
          nested_threads=(MAX(parai%ncpus-1,1))
#if !defined(_INTEL_MKL)
          CALL stopgm(procedureN, 'Overlapping communication and computation: Behavior of &
               BLAS routine inside parallel region not checked',&
               __LINE__,__FILE__)
#endif
       ELSE
          nthreads=1
          nested_threads=parai%ncpus
       END IF
       methread=0
       IF(buffcount.EQ.1)THEN
          nthreads=1
          nested_threads=parai%ncpus
       END IF
       methread=0
       IF(buffcount.eq.1)THEN
          nthreads=1
          nested_threads=parai%ncpus
       END IF
       methread=0
       !$ locks=.TRUE.
       !$ locks(1)=.FALSE.
       !$omp parallel if(nthreads.eq.2) num_threads(nthreads) &
       !$omp private(methread,buff,ld_dai,end_dai,start_dai)&
       !$omp shared(locks,nthreads) proc_bind(close)
       !$ methread = omp_get_thread_num()
       !$ IF (nested_threads .NE. parai%ncpus) THEN
       !$    IF (methread .EQ. 1 .OR. nthreads .EQ. 1) THEN
       !$       CALL omp_set_num_threads(nested_threads)
#ifdef _INTEL_MKL
       !$       CALL mkl_set_dynamic(0)
#endif
       !$       CALL omp_set_nested(.TRUE.)
       !$    END IF
       !$ END IF
       !$OMP barrier

       IF(methread.EQ.1.OR.nthreads.EQ.1) THEN
          DO buff=2,buffcount
             ld_dai=ld_buffer(buff)
             start_dai=start_buffer(buff)
             end_dai=start_dai-1+ld_dai*nstate
             IF(ld_dai.GT.0)THEN
                CALL proj_beta(na_buff(:,:,buff),igeq0,nstate,c0,nkpt%ngwk,eigkr(:,:,ikind),&
                     twnl(:,:,:,ikind),eiscr,t,nkpt%ngwk,1,dai(start_dai:end_dai),&
                     ld_dai/imagp,.TRUE.,tkpts%tkpnt,geq0,gktemp)
             END IF
             !$ locks(buff) = .FALSE.
             !$omp flush(locks)
          ENDDO
       ENDIF

       IF(methread.EQ.0.or.nthreads.EQ.1) THEN
          DO buff=1,buffcount
             ld_dai=ld_buffer(buff)
             start_dai=start_buffer(buff)
             end_dai=start_dai-1+ld_dai*nstate
             CALL TISET(procedureN//'_barrier',ISUB2)
             !$omp flush(locks)
             !$ DO WHILE ( locks(buff) )
             !$omp flush(locks)
             !$ END DO
             CALL TIHALT(procedureN//'_barrier',ISUB2)
             IF(autotune_it.GT.0.AND.autotune_it.LE.cnti%rnlsm_autotune_maxit) temp=m_walltime()
             IF(ld_dai.GT.0)THEN
                CALL mp_sum(dai(start_dai:end_dai),end_dai-start_dai+1,parai%allgrp)
             END IF
             IF(autotune_it.GT.0.AND.autotune_it.LE.cnti%rnlsm_autotune_maxit)&
                  timings(2)=timings(2)+m_walltime()-temp
          ENDDO
       ENDIF

       !$omp barrier
       !$ IF (nested_threads .NE. parai%ncpus) THEN
       !$    IF (methread .EQ. 1 .OR. nthreads .EQ. 1) THEN
       !$       CALL omp_set_num_threads(parai%ncpus)
#ifdef _INTEL_MKL
       !$       CALL mkl_set_dynamic(1)
#endif
       !$       CALL omp_set_nested(.FALSE.)
       !$    END IF
       !$ END IF

       !$OMP end parallel
       IF(buffcount.GT.1)THEN
          CALL sort_dfnl(buffcount,na_buff(:,:,:),dai,dfnl_packed,&
               start_buffer,ld_buffer)
       END IF
       !only local states need to be unpacked
       IF(tkpts%tkpnt)THEN
          CALL unpack_dfnl(na_grp(:,:,parai%cp_inter_me),&
               dfnl_packed(:,parap%nst12(parai%mepos,1):parap%nst12(parai%mepos,2)),&
               unpacked_k=dfnl(:,:,:,:,:,ikind))
       ELSE
          CALL unpack_dfnl(na_grp(:,:,parai%cp_inter_me),&
               dfnl_packed(:,parap%nst12(parai%mepos,1):parap%nst12(parai%mepos,2)),&
               unpacked=dfnla)
       END IF
       IF(parai%cp_nogrp.GT.1)CALL cp_grp_redist_dfnl_fnl(.FALSE.,.TRUE.,nstate,ikind)
       !set buffcount/buffsizes
       IF(cnti%rnlsm_autotune_maxit.GT.0)THEN
          IF(autotune_it.EQ.cnti%rnlsm_autotune_maxit)THEN
             !stop timemeasurements
             autotune_it=autotune_it+1
             IF(paral%parent)THEN
                CALL tune_rnlsm(cnti%rnlsm2_bc,cntr%rnlsm2_b1,&
                     cntr%rnlsm2_b2,timings,INT(il_dfnl_packed(1),KIND=int_4),nstate)
                WRITE(6,*) "####rnlsm2_autotuning results####"
                WRITE(6,*) cnti%rnlsm2_bc
                WRITE(6,*) cntr%rnlsm2_b1
                WRITE(6,*) cntr%rnlsm2_b2
                WRITE(6,*) timings
             END IF
             !broadcast results
             CALL mp_bcast(cnti%rnlsm2_bc,parai%source,parai%allgrp)
             CALL mp_bcast(cntr%rnlsm2_b1,parai%source,parai%allgrp)
             CALL mp_bcast(cntr%rnlsm2_b2,parai%source,parai%allgrp)
          END IF
       END IF
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_t,t,procedureN//'_t')
    CALL free_scratch(il_gktemp,gktemp,procedureN//'_gktemp')
    CALL free_scratch(il_eiscr,eiscr,procedureN//'_eiscr')
#else
    DEALLOCATE(eiscr, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate eiscr',&
         __LINE__,__FILE__)
    DEALLOCATE(gktemp, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate gktemp',&
         __LINE__,__FILE__)
    DEALLOCATE(t, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate t',&
         __LINE__,__FILE__)
#endif
    IF(buffcount.GT.1)THEN
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_dai,dai,procedureN//'_dai')
#else
       DEALLOCATE(dai, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate dai',&
            __LINE__,__FILE__)
#endif
    END IF
    DEALLOCATE(na_buff, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate na_buff',&
         __LINE__,__FILE__)
    DEALLOCATE(na_grp, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate na_grp',&
         __LINE__,__FILE__)

    CALL TIHALT(procedureN,ISUB)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE RNLSM2
  ! ==================================================================

END MODULE rnlsm2_utils
