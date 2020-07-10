MODULE csize_utils
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes
  USE dotp_utils,                      ONLY: dotp_c1_cp
  USE elct,                            ONLY: crge
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_max,&
                                             mp_sum
  USE parac,                           ONLY: parai
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             nkpt,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
#ifdef __PARALLEL
  USE mpi_f08
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: csize

CONTAINS

  ! ==================================================================
  SUBROUTINE csize(c2,nstate,gemax,cnorm,use_cp_grps,special)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate
    COMPLEX(real_8),INTENT(IN)               :: c2(nkpt%ngwk,nstate)
    REAL(real_8),INTENT(OUT)                 :: gemax, cnorm
    LOGICAL,INTENT(IN),OPTIONAL              :: use_cp_grps, special

    CHARACTER(*), PARAMETER                  :: procedureN = 'csize'

#ifdef __PARALLEL
    INTEGER                                  :: i, iiabs, nocc, isub, &
                                                ngwk_local, ibeg_c0
    type(MPI_COMM)                           :: gid
#else
    INTEGER                                  :: i, iiabs, nocc, isub, &
                                                ngwk_local, ibeg_c0, gid
#endif

    LOGICAL                                  :: geq0_local, cp_active, sp
    REAL(real_8), EXTERNAL                   :: ddot

    INTEGER, EXTERNAL                        :: izamax
    ! ==--------------------------------------------------------------==
    !TK in noforce cnorm and gemax were multiplied by the occupation number
    ! special case to get it done in the same way...
    CALL tiset(procedureN,isub)
    IF(PRESENT(special))THEN
       sp=special
    ELSE
       sp=.FALSE.
    END IF
    gemax=0.0_real_8
    cnorm=0.0_real_8
    nocc=0.0_real_8
    IF(PRESENT(use_cp_grps))THEN
       cp_active=use_cp_grps
    ELSE
       cp_active=.FALSE.
    END IF
    IF(cp_active)THEN
       CALL cp_grp_get_sizes(ngwk_l=ngwk_local,firstk_g=ibeg_c0,geq0_l=geq0_local)
       gid=parai%cp_grp
    ELSE
       ngwk_local=nkpt%ngwk
       ibeg_c0=1
       geq0_local=geq0
       gid=parai%allgrp
    END IF
    IF(ngwk_local.GT.0)THEN
       IF (tkpts%tkpnt) THEN
          !$omp parallel do private(i,iiabs) reduction(max:gemax) reduction(+:cnorm,nocc)
          DO i=1,nstate
             IF(crge%f(i,1).LE.1.e-5_real_8) CYCLE
             nocc=nocc+1
             IF(sp)THEN
                cnorm=cnorm+crge%f(i,1)*ddot(2*ngwk_local,c2(ibeg_c0,i),1,c2(ibeg_c0,i),1)
                iiabs=izamax(ngwk_local,c2(ibeg_c0,i),1)
                gemax=MAX(ABS(crge%f(i,1)*c2(iiabs,i)),gemax)
             ELSE
                cnorm=cnorm+ddot(2*ngwk_local,c2(ibeg_c0,i),1,c2(ibeg_c0,i),1)
                iiabs=izamax(ngwk_local,c2(ibeg_c0,i),1)
                gemax=MAX(ABS(c2(iiabs,i)),gemax)
             END IF
          END DO
       ELSE
          !$omp parallel do private(i,iiabs) reduction(max:gemax) reduction(+:cnorm,nocc)
          DO i=1,nstate
             IF(crge%f(i,1).LE.1.e-5_real_8) CYCLE
             nocc=nocc+1
             IF(sp)THEN
                cnorm=cnorm+crge%f(i,1)*dotp_c1_cp(ngwk_local,c2(ibeg_c0,i),geq0_local)
                iiabs=izamax(ngwk_local,c2(ibeg_c0,i),1)
                gemax=MAX(ABS(crge%f(i,1)*c2(iiabs,i)),gemax)
             ELSE
                cnorm=cnorm+dotp_c1_cp(ngwk_local,c2(ibeg_c0,i),geq0_local)
                iiabs=izamax(ngwk_local,c2(ibeg_c0,i),1)
                gemax=MAX(ABS(c2(iiabs,i)),gemax)
             END IF
          END DO
       END IF
    END IF
    CALL mp_sum(cnorm,gid)
    CALL mp_max(gemax,gid)
    cnorm=SQRT(cnorm/REAL(nocc*spar%ngwks,kind=real_8))
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE csize
  ! ==================================================================

END MODULE csize_utils
