MODULE beta_utils
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nlps_com,&
                                             nghtol
  USE system,                          ONLY: maxsys

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: build_beta
  PUBLIC :: build_beta_deriv

CONTAINS
  ! Author: Tobias Kloeffel, FAU Erlangen-Nuernberg April 2019
  ! Contains routines to build betaprojector arrays as needed in nlforce,
  ! rnlsm1/2 spsi
  ! WARNING: Call from omp parallel region to get omp parallelism here
  ! WARNING: no omp barrier at the end of both routines
  ! WARNING: geq0 component not considered here
  SUBROUTINE build_beta_deriv(na,eigkr,gktemp,twnl,eiscr,t,ngw,ig_start,ld_eiscr)
    ! ==--------------------------------------------------------------==
    ! == build derivative of beta projectors according to na          ==
    ! == call from parallel region to get omp parallelism here        ==
    ! == WARNING: do not forget to set geq0 component outside!        ==
    ! == WARNING: no omp barrier at the end of this routine!          ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: na(2,*),ngw,ig_start,ld_eiscr
    COMPLEX(real_8),INTENT(IN)               :: eigkr(ngw,*)
    REAL(real_8),INTENT(IN)                  :: gktemp(ld_eiscr,*), &
                                                twnl(ngw,maxsys%nhxs,ions1%nsp,*)
    COMPLEX(real_8),INTENT(OUT)              :: eiscr(ld_eiscr,*)
    REAL(real_8),INTENT(OUT)                 :: t(ld_eiscr,*)

    INTEGER                                  :: isa0, is, start_ia, end_ia, ia_sum, &
                                                iv, k, ia, isa, offset, offset0, ig
    REAL(real_8)                             :: cir,cii
    COMPLEX(real_8)                          :: ci
    ! ==--------------------------------------------------------------==
    offset0=0
    isa0=0
    DO is=1,ions1%nsp
       start_ia=na(1,is)
       end_ia=na(2,is)
       ia_sum=end_ia-start_ia+1
       IF(ia_sum.GT.0)THEN
          DO iv=1,nlps_com%ngh(is)
             ci=(0.0_real_8,-1.0_real_8)**(nghtol(iv,is)+1)
             cir=REAL(ci,kind=real_8)
             cii=AIMAG(ci)
             !Make use of the special structure of CI
             IF (ABS(cir).GT.0.5_real_8) THEN
                !CI is real
                DO k=1,3
                   !$omp barrier
                   !$omp do
                   DO ig=1,ld_eiscr
                      t(ig,1)=twnl(ig+ig_start-1,iv,is,1)*cir*gktemp(ig,k)
                   END DO
                   !$OMP do
                   DO ia=start_ia,end_ia
                      isa=isa0+ia
                      offset=offset0+ia-start_ia+1
                      CALL build_betaproj_r(eiscr(:,offset),eigkr(ig_start,isa),t,ld_eiscr)
                   END DO
                   !$omp end do nowait
                   offset0=offset0+ia_sum
                END DO
             ELSE
                !CI is imaginary
                DO k=1,3
                   !$omp barrier
                   !$omp do
                   DO ig=1,ld_eiscr
                      t(ig,1)=twnl(ig+ig_start-1,iv,is,1)*cii*gktemp(ig,k)
                   END DO
                   !$OMP do
                   DO ia=start_ia,end_ia
                      isa=isa0+ia
                      offset=offset0+ia-start_ia+1
                      CALL build_betaproj_i(eiscr(:,offset),eigkr(ig_start,isa),t,ld_eiscr)
                   END DO
                   !$omp end do nowait
                   offset0=offset0+ia_sum
                END DO
             END IF
          END DO
       END IF
       isa0=isa0+ions0%na(is)
    END DO
  END SUBROUTINE build_beta_deriv
  ! ==================================================================
  SUBROUTINE build_beta(na,eigkr,twnl,eiscr,t,ngw,ig_start,ld_eiscr)
    ! ==--------------------------------------------------------------==
    ! == build derivative of beta projectors according to na          ==
    ! == call from parallel region to get omp parallelism here        ==
    ! == WARNING: do not forget to set geq0 component outside!        ==
    ! == WARNING: no omp barrier at the end of this routine!          ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: na(2,*),ngw,ig_start,ld_eiscr
    COMPLEX(real_8),INTENT(IN)               :: eigkr(ngw,*)
    REAL(real_8),INTENT(IN)                  :: twnl(ngw,maxsys%nhxs,ions1%nsp,*)
    COMPLEX(real_8),INTENT(OUT)              :: eiscr(ld_eiscr,*)
    REAL(real_8),INTENT(OUT)                 :: t(ld_eiscr,*)
    INTEGER                                  :: isa0, is, start_ia, end_ia, ia_sum, iv,&
                                                ia, isa, offset, offset0, ig
    REAL(real_8)                             :: cir,cii
    COMPLEX(real_8)                          :: ci
    ! ==--------------------------------------------------------------==
    offset0=0
    isa0=0
    DO is=1,ions1%nsp
       start_ia=na(1,is)
       end_ia=na(2,is)
       ia_sum=end_ia-start_ia+1
       IF(ia_sum.GT.0)THEN
          DO iv=1,nlps_com%ngh(is)
             ci=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
             cir=REAL(ci,kind=real_8)
             cii=AIMAG(ci)
             !$omp barrier
             !Make use of the special structure of CI
             IF (ABS(cir).GT.0.5_real_8) THEN
                !CI is real
                !$omp do
                DO ig=1,ld_eiscr
                   t(ig,1)=twnl(ig+ig_start-1,iv,is,1)*cir
                END DO
                !$OMP do
                DO ia=start_ia,end_ia
                   isa=isa0+ia
                   offset=offset0+ia-start_ia+1
                   CALL build_betaproj_r(eiscr(:,offset),eigkr(ig_start,isa),t,ld_eiscr)
                END DO
                !$omp end do nowait
             ELSE
                !CI is imaginary
                !$omp do
                DO ig=1,ld_eiscr
                   t(ig,1)=twnl(ig+ig_start-1,iv,is,1)*cii
                END DO
                !$OMP do
                DO ia=start_ia,end_ia
                   isa=isa0+ia
                   offset=offset0+ia-start_ia+1
                   CALL build_betaproj_i(eiscr(:,offset),eigkr(ig_start,isa),t,ld_eiscr)
                END DO
                !$omp end do nowait
             END IF
             offset0=offset0+ia_sum
          END DO
       END IF
       isa0=isa0+ions0%na(is)
    END DO
    RETURN
  END SUBROUTINE build_beta
  ! ==================================================================
  PURE SUBROUTINE build_betaproj_r(eiscr,eigr_,twnl_,lda)
    INTEGER,INTENT(IN)                       :: lda
    COMPLEX(real_8), INTENT(OUT)             :: eiscr(lda)
    COMPLEX(real_8),INTENT(IN)               :: eigr_(lda,*)
    REAL(real_8),INTENT(IN)                  :: twnl_(lda,*)
    REAL(real_8)                             :: t1, t2
    INTEGER                                  :: ig

    DO ig=1,lda
       t1=REAL(eigr_(ig,1))*twnl_(ig,1)
       t2=AIMAG(eigr_(ig,1))*twnl_(ig,1)
       eiscr(ig)=CMPLX(t1,t2,kind=real_8)
    END DO

  END SUBROUTINE build_betaproj_r
  ! ==================================================================
  PURE SUBROUTINE build_betaproj_i(eiscr,eigr_,twnl_,lda)
    INTEGER,INTENT(IN)                       :: lda
    COMPLEX(real_8),INTENT(OUT)              :: eiscr(lda)
    COMPLEX(real_8),INTENT(IN)               :: eigr_(lda,*)
    REAL(real_8),INTENT(IN)                  :: twnl_(lda,*)
    REAL(real_8)                             :: t1, t2
    INTEGER                                  :: ig

    DO ig=1,lda
       t1=REAL(eigr_(ig,1))*twnl_(ig,1)
       t2=AIMAG(eigr_(ig,1))*twnl_(ig,1)
       eiscr(ig)=CMPLX(-t2,t1,kind=real_8)
    END DO

  END SUBROUTINE build_betaproj_i
  ! ==================================================================
END MODULE beta_utils
