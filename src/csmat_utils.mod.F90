#include "cpmd_global.h"

MODULE csmat_utils
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms
  USE cvan,                            ONLY: qq
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE distribution_utils,              ONLY: dist_atoms
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: nlps_com
  USE nort,                            ONLY: nort_com,&
                                             nort_ovlap
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: fnl_packed
  USE spin,                            ONLY: spin_mod
  USE summat_utils,                    ONLY: summat
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             maxsys,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: csmat

CONTAINS

  ! ==================================================================
  ! == FOR TKPNT=.TRUE. FNL IS COMPLEX -- CSMAT_C WILL BE WRITTEN   ==
  ! ==================================================================
  SUBROUTINE csmat(a,c0,nstate,ikind,full,store_nonort,only_parent)
    ! ==--------------------------------------------------------------==
    ! ==         COMPUTES THE OVERLAP MATRIX A = < C0 |S| C0 >        ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: ikind, nstate
    COMPLEX(real_8),INTENT(IN)               :: c0(ncpw%ngw,nstate)
    REAL(real_8),INTENT(OUT)                 :: a(nstate,*)
    LOGICAL,INTENT(IN)                       :: full,only_parent,store_nonort
    INTEGER                                  :: i, j, ia, is, ispin, nspin,&
                                                ia_sum, tot_work, ns(2), nmin(2), &
                                                off_i, off_mat, off_fnl, ia_fnl, &
                                                start_fnl, end_fnl, fnl_start, start_mat, &
                                                end_mat, il_fnlat(3), il_fnlatj(3), &
                                                na(2,ions1%nsp),na_fnl(2,ions1%nsp),&
                                                na_grp(2,ions1%nsp,0:parai%cp_nogrp-1),&
                                                isub, ierr
    REAL(real_8)                             :: fractions(parai%nproc),selem, temp
    CHARACTER(*), PARAMETER                  :: procedureN = 'csmat'
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8),POINTER __CONTIGUOUS        :: fnlat(:,:,:), fnlatj(:,:,:)
#else
    REAL(real_8),ALLOCATABLE                 :: fnlat(:,:,:), fnlatj(:,:,:)
#endif
    INTEGER,ALLOCATABLE,SAVE                 :: na_buff(:,:,:)

    CALL tiset(procedureN,isub)
    IF (cntl%tfdist) CALL stopgm('CSMAT','TFDIST NOT IMPLEMENTED',&
         __LINE__,__FILE__)

    !we only need the upper part here
    CALL ovlap(nstate,a,c0,c0,redist=.FALSE.,full=.FALSE.)
    IF(store_nonort) CALL store_ovlap(a,nstate)

    IF (pslo_com%tivan) THEN
       CALL cp_grp_split_atoms(na_grp)
       na_fnl(:,:)=na_grp(:,:,parai%cp_inter_me)
       IF (.NOT.ALLOCATED(na_buff))THEN
          !distribute atoms between procs
          ALLOCATE(na_buff(2,ions1%nsp,0:parai%nproc-1), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate na_buff',&
               __LINE__,__FILE__)
          fractions=1.0_real_8/REAL(parai%nproc,KIND=real_8)
          CALL dist_atoms(parai%nproc,fractions,na_fnl,na_buff,only_uspp=.TRUE.)
       END IF
       na(:,:)=na_buff(:,:,parai%mepos)

       ! VANDERBILT PP
       ns(1)=nstate
       ns(2)=0
       nmin(1)=1
       nmin(2)=0
       nspin=1
       IF (cntl%tlsd) THEN
          ns(1)=spin_mod%nsup
          ns(2)=spin_mod%nsdown
          nmin(1)=1
          nmin(2)=spin_mod%nsup+1
          nspin=2
       ENDIF
       !end preparation

       !calculate local part of total work
       tot_work=0
       DO is=1,ions1%nsp
          tot_work=tot_work+(na(2,is)-na(1,is)+1)*nlps_com%ngh(is)
       END DO

       il_fnlat(1)=tot_work
       il_fnlat(2)=MAXVAL(ns)
       il_fnlat(3)=nspin

       il_fnlatj(1)=tot_work
       il_fnlatj(2)=MAXVAL(ns)
       il_fnlatj(3)=nspin
       IF(tot_work.GT.0)THEN
#ifdef _USE_SCRATCHLIBRARY
          CALL request_scratch(il_fnlat,fnlat,procedureN//'_fnlat')
          CALL request_scratch(il_fnlatj,fnlatj,procedureN//'_fnlatj')
#else
          ALLOCATE(fnlat(il_fnlat(1),il_fnlat(2),il_fnlat(3)), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlat',&
               __LINE__,__FILE__)
          ALLOCATE(fnlatj(il_fnlatj(1),il_fnlatj(2),il_fnlatj(3)), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlatj',&
               __LINE__,__FILE__)
#endif
          !$omp parallel private(ispin,off_i,i,off_mat,off_fnl,is,ia_fnl,ia_sum,&
          !$omp start_fnl,end_fnl,fnl_start,start_mat,end_mat)
          DO ispin=1,nspin
             off_i=0
             IF(ispin.EQ.2) off_i=spin_mod%nsup
             !$omp do schedule(static)
             DO i=1,ns(ispin)
                off_mat=0
                off_fnl=0
                DO is=1,ions1%nsp
                   ia_fnl=na_fnl(2,is)-na_fnl(1,is)+1
                   ia_sum=na(2,is)-na(1,is)+1
                   !nc alreday filtered out in dist_atoms
                   IF(ia_sum.GT.0)THEN
                      !starting index fnl_packed
                      fnl_start=na(1,is)-na_fnl(1,is)
                      !fnl_range
                      start_fnl=off_fnl+1
                      end_fnl=off_fnl+ia_fnl*nlps_com%ngh(is)
                      !mat_range
                      start_mat=off_mat+1
                      end_mat=off_mat+ia_sum*nlps_com%ngh(is)
                      CALL prepare_matrix(fnl_packed(start_fnl:end_fnl,i+off_i),&
                           fnlat(start_mat:end_mat,i,ispin),&
                           fnlatj(start_mat:end_mat,i,ispin),qq(:,:,is),nlps_com%ngh(is),&
                           ia_sum,ia_fnl,fnl_start)
                      off_mat=off_mat+ia_sum*nlps_com%ngh(is)
                   END IF
                   off_fnl=off_fnl+ia_fnl*nlps_com%ngh(is)
                END DO
             END DO
          END DO
          !$omp end parallel
          DO ispin=1,nspin
#ifdef _HAS_DGEMMT
             CALL dgemmt('U','T','N',ns(ispin),tot_work,1.0_real_8,&
                  fnlat(1,1,ispin),tot_work,fnlatj(1,1,ispin),tot_work,1.0_real_8,&
                  a(nmin(ispin),nmin(ispin)),nstate)
#else
             CALL dgemm('T','N',ns(ispin),ns(ispin),tot_work,1.0_real_8,&
                  fnlat(1,1,ispin),tot_work,fnlatj(1,1,ispin),tot_work,1.0_real_8,&
                  a(nmin(ispin),nmin(ispin)),nstate)
#endif
          END DO
#ifdef _USE_SCRATCHLIBRARY
          CALL free_scratch(il_fnlatj,fnlatj,procedureN//'_fnlatj')
          CALL free_scratch(il_fnlat,fnlat,procedureN//'_fnlat')
#else
          DEALLOCATE(fnlat, stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fnlat',&
               __LINE__,__FILE__)
          DEALLOCATE(fnlatj, stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fnlatj',&
               __LINE__,__FILE__)
#endif
       END IF
    END IF
    CALL summat(a,nstate,symmetrization=full,lsd=.TRUE.,gid=parai%cp_grp,&
         parent=only_parent)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE csmat
  ! ==================================================================
  SUBROUTINE store_ovlap(a,nstate)
    ! ==--------------------------------------------------------------==
    ! ==         STORE OVLAP MATRIX FOR LATER USE IN CROTWF           ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate
    REAL(real_8),INTENT(IN)                  :: a(nstate,*)
    INTEGER                                  :: i,j, ierr
    REAL(real_8)                             :: temp, selem
    CHARACTER(*), PARAMETER                  :: procedureN = 'store_nort'
    ! ==--------------------------------------------------------------==
    IF (.NOT. ALLOCATED(nort_ovlap)) THEN
       ALLOCATE(nort_ovlap(nstate,nstate), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate nort_ovlap',&
            __LINE__,__FILE__)
    ELSE
       IF(SIZE(nort_ovlap).NE.nstate*nstate)THEN
          DEALLOCATE(nort_ovlap, stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate nort_ovlap',&
               __LINE__,__FILE__)
          ALLOCATE(nort_ovlap(nstate,nstate), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate nort_ovlap',&
               __LINE__,__FILE__)
       END IF
    END IF
    IF(cntl%tlsd)THEN
       !$omp parallel private(i,j)
       !$omp do
       DO i=1,spin_mod%nsup
          DO j=1,i
             nort_ovlap(j,i)=a(j,i)
          END DO
       END DO
       !$omp end do nowait
       !$omp  do
       DO i=spin_mod%nsup+1,nstate
          DO j=spin_mod%nsup+1,i
             nort_ovlap(j,i)=a(j,i)
          END DO
       END DO
       !$omp end parallel
    ELSE
       !$omp parallel do private(i,j)
       DO i=1,nstate
          DO j=1,i
             nort_ovlap(j,i)=a(j,i)
          END DO
       END DO
    END IF
    !lapack only needs upper or lower part
    CALL summat(nort_ovlap,nstate,symmetrization=.FALSE.,lsd=.TRUE.,gid=parai%cp_grp,&
         parent=.FALSE.)
    !always check threshold.
    temp=0.0_real_8
    IF(cntl%tlsd)THEN
       !$omp parallel private(i,j,selem)reduction(max:temp)
       !$omp do
       DO i=1,spin_mod%nsup
          DO j=1,i-1
             selem=ABS(nort_ovlap(j,i))
             IF (temp.LT.selem) temp=selem
          ENDDO
       ENDDO
       !$omp end do nowait
       !$omp do
       DO i=spin_mod%nsup+1,nstate
          DO j=spin_mod%nsup+1,i-1
             selem=ABS(nort_ovlap(j,i))
             IF (temp.LT.selem) temp=selem
          ENDDO
       ENDDO
       !$omp end parallel
    ELSE
       !$omp parallel do private(i,j,selem)reduction(max:temp)
       DO i=1,nstate
          DO j=1,i-1
             selem=ABS(nort_ovlap(j,i))
             IF (temp.LT.selem) temp=selem
          ENDDO
       ENDDO
    END IF
    nort_com%scond=temp

    RETURN
  END SUBROUTINE store_ovlap
  ! ==================================================================
  PURE SUBROUTINE prepare_matrix(fnl_p,fnli,fnlj,qq_,ngh,ia_sum,ia_fnl,fnl_start)
    INTEGER,INTENT(IN)                       :: ngh,ia_sum,ia_fnl,fnl_start
    REAL(real_8),INTENT(IN)                  :: fnl_p(ia_fnl,ngh,*)
    REAL(real_8),INTENT(IN) __CONTIGUOUS     :: qq_(:,:)
    REAL(real_8),INTENT(OUT)                 :: fnli(ia_sum,ngh,*),fnlj(ia_sum,ngh,*)
    INTEGER                                  :: iv,ia,jv

    DO iv=1,ngh
       DO ia=1,ia_sum
          fnli(ia,iv,1)=fnl_p(ia+fnl_start,iv,1)
          fnlj(ia,iv,1)=0.0_real_8
       END DO
    END DO
    DO iv=1,ngh
       DO jv=1,ngh
          IF (ABS(qq_(jv,iv)).GT.1.e-5_real_8) THEN
             DO ia=1,ia_sum
                fnlj(ia,iv,1)=fnlj(ia,iv,1)&
                     +qq_(jv,iv)*fnli(ia,jv,1)
             END DO
          END IF
       END DO
    END DO
  END SUBROUTINE prepare_matrix
  ! ==================================================================
END MODULE csmat_utils
