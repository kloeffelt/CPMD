MODULE crotwf_utils
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai
  USE rotate_utils,                    ONLY: rotate
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE utils,                           ONLY: dspevy
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: crotwf
  PUBLIC :: give_scr_crotwf

CONTAINS

  ! ==================================================================
  SUBROUTINE crotwf(c0,cm,c2,sc0,nstate,gam,use_cp_grps)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate
    COMPLEX(real_8),INTENT(INOUT)            :: sc0(ncpw%ngw,nstate), &
                                                c2(ncpw%ngw,nstate), &
                                                cm(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    REAL(real_8),INTENT(OUT)                 :: gam(nstate,nstate)
    LOGICAL, INTENT(IN)                      :: use_cp_grps

    CHARACTER(*), PARAMETER                  :: procedureN = 'crotwf'

    INTEGER                                  :: i, ic1b, ierr, iopt, j, k,&
                                                ibeg, iend, ig
    REAL(real_8), ALLOCATABLE                :: aux(:), c1(:), w(:)

    ALLOCATE(c1(nstate*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(c1)
    ALLOCATE(w(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(aux(3*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL ovlap(nstate,gam,c0,c0)
    IF (.NOT.cntl%tlsd) THEN
       k=1
       DO j=1,nstate
          DO i=j,nstate
             c1(k)=gam(i,j)
             k=k+1
          ENDDO
       ENDDO
       CALL mp_sum(c1,nstate*nstate,parai%allgrp)
       iopt=1
       CALL dspevy(iopt,c1,w,gam,nstate,nstate,aux,3*nstate)
    ELSE
       CALL mp_sum(gam,nstate*nstate,parai%allgrp)
       k=1
       DO j=1,spin_mod%nsup
          DO i=j,spin_mod%nsup
             c1(k)=gam(i,j)
             k=k+1
          ENDDO
       ENDDO
       iopt=1
       CALL dspevy(iopt,c1,w,gam,spin_mod%nsup,spin_mod%nsup,aux,3*spin_mod%nsup)
       CALL dcopy(spin_mod%nsup*spin_mod%nsup,gam(1,1),1,c1,1)
       ic1b = spin_mod%nsup*spin_mod%nsup
       k=1
       DO j=1,spin_mod%nsdown
          DO i=j,spin_mod%nsdown
             c1(ic1b+k)=gam(spin_mod%nsup+i,spin_mod%nsup+j)
             k=k+1
          ENDDO
       ENDDO
       iopt=1
       CALL dspevy(iopt,c1(ic1b+1),w,gam,spin_mod%nsdown,spin_mod%nsdown,aux,3*spin_mod%nsdown)
       CALL dcopy(spin_mod%nsdown*spin_mod%nsdown,gam(1,1),1,c1(ic1b+1),1)
       CALL zeroing(gam)!,nstate*nstate)
       k=1
       DO j=1,spin_mod%nsup
          DO i=1,spin_mod%nsup
             gam(i,j)=c1(k)
             k=k+1
          ENDDO
       ENDDO
       k=1
       DO j=1,spin_mod%nsdown
          DO i=1,spin_mod%nsdown
             gam(spin_mod%nsup+i,spin_mod%nsup+j)=c1(ic1b+k)
             k=k+1
          ENDDO
       ENDDO
    ENDIF

    ! to avoid problems we bcast the result (we dont need the eigvals)
    CALL mp_bcast(gam,nstate**2,parai%io_source,parai%cp_grp)
    IF(use_cp_grps)THEN
       CALL cp_grp_get_sizes(first_g=ibeg,last_g=iend)
    ELSE
       ibeg=1
       iend=ncpw%ngw
    END IF

    CALL rotate(1.0_real_8,c0,0.0_real_8,sc0,gam,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown,&
         redist=.NOT.use_cp_grps,use_cp=use_cp_grps)
    !$omp parallel do private(ig,i)
    DO i=1,nstate
       DO ig=ibeg,iend
          c0(ig,i)=sc0(ig,i)
       END DO
    END DO
    CALL rotate(1.0_real_8,cm,0.0_real_8,sc0,gam,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown,&
         redist=.NOT.use_cp_grps,use_cp=use_cp_grps)
    !$omp parallel do private(ig,i)
    DO i=1,nstate
       DO ig=ibeg,iend
          cm(ig,i)=sc0(ig,i)
       END DO
    END DO
    CALL rotate(1.0_real_8,c2,0.0_real_8,sc0,gam,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown,&
         redist=.NOT.use_cp_grps,use_cp=use_cp_grps)
    !$omp parallel do private(ig,i)
    DO i=1,nstate
       DO ig=ibeg,iend
          c2(ig,i)=sc0(ig,i)
       END DO
    END DO
    ! ==--------------------------------------------------------------==
    DEALLOCATE(c1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(w,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE crotwf
  ! ==================================================================
  SUBROUTINE give_scr_crotwf(lcrotwf,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lcrotwf
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

! ==--------------------------------------------------------------==

    lcrotwf=nstate*nstate+4*nstate
    tag   ='NSTATE*NSTATE+4*NSTATE'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_crotwf
  ! ==================================================================

END MODULE crotwf_utils
