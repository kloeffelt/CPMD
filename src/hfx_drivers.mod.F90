MODULE hfx_drivers
  USE hfx_utils,                       ONLY: hfx_old,&
                                             hfxpsi_old,&
                                             hfxrpa_old, &
                                             hfx_scdm, hfx_ace !SM
  USE hfxmod,                          ONLY: hfxc3
  USE kinds,                           ONLY: real_8
  USE pw_hfx,                          ONLY: hfx_new
  USE pw_hfx_resp,                     ONLY: hfxpsi_new,&
                                             hfxrpa_incore,&
                                             hfxrpa_new
  USE pw_hfx_resp_types,               ONLY: hfx_resp
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE ace_hfx  !SM

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hfx
  PUBLIC :: hfxpsi
  PUBLIC :: hfxrpa


CONTAINS


  ! ==================================================================
  SUBROUTINE hfxrpa(c0,c1,c2,psia,nstate,tcis)
    ! ==================================================================
    ! == HARTREE-FOCK EXCHANGE CONTRIBUTION TO RPA                    ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c1(:,:), c2(:,:), &
                                                psia(:)
    INTEGER                                  :: nstate
    LOGICAL                                  :: tcis

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxrpa'

    INTEGER                                  :: isub

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF (hfxc3%use_new_hfx) THEN
       IF (hfx_resp%is_set.AND.hfx_resp%init) THEN
          CALL hfxrpa_incore(c0,c1,c2,psia,nstate,tcis)
       ELSE
          CALL hfxrpa_new(c0,c1,c2,psia,nstate,tcis)
       ENDIF
    ELSE
       CALL hfxrpa_old(c0,c1,c2,psia,nstate,tcis)
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxrpa
  ! ==================================================================
  SUBROUTINE hfxpsi(c0,cpsi,c2,f,sign,psia,nstate,norb)
    ! ==================================================================
    ! == HARTREE-FOCK EXCHANGE                                        ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), cpsi(:,:), c2(:,:)
    REAL(real_8)                             :: f(:), sign
    COMPLEX(real_8)                          :: psia(:)
    INTEGER                                  :: nstate, norb

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxpsi'

    INTEGER                                  :: isub

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF (hfxc3%use_new_hfx) THEN
       CALL hfxpsi_new(c0,cpsi,c2,f,sign,psia,nstate,norb)
    ELSE
       CALL hfxpsi_old(c0,cpsi,c2,f,sign,psia,nstate,norb)
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxpsi
  ! ==================================================================
  SUBROUTINE hfx(c0,c2,f,psia,nstate,ehfx,vhfx,redist_c2)
    ! ==================================================================
    ! == HARTREE-FOCK EXCHANGE WRAPPER FOR OLD/NEW MODULES            ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:)
    REAL(real_8)                             :: f(:)
    COMPLEX(real_8)                          :: psia(:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: ehfx, vhfx
    LOGICAL                                  :: redist_c2

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfx'

    INTEGER                                  :: isub

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF (hfxc3%use_new_hfx) THEN
     !
     if(use_ace.and.status_ace)then
       !write(6,*)"redist_c2=",redist_c2 
       if(switch_ace)then  !SM
         call hfx_ace(c0,c2,f,psia,nstate,ehfx,vhfx)
       elseif(hfx_scdm_status)then
          CALL hfx_scdm(c0,c2,f,psia,nstate,ehfx,vhfx,redist_c2) !TODO SM
       else
           CALL hfx_new(c0,c2,f,psia,nstate,ehfx,vhfx,redist_c2)
       endif
     else 
       CALL hfx_new(c0,c2,f,psia,nstate,ehfx,vhfx,redist_c2)
     endif
     !
    ELSE
       CALL hfx_old(c0,c2,f,psia,nstate,ehfx,vhfx)
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfx
  ! ==================================================================

END MODULE hfx_drivers
