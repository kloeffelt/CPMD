#include "cpmd_global.h"
  
SUBROUTINE cpmd_dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  USE kinds,                          ONLY: real_8
  CHARACTER(1), INTENT(IN)    :: TRANSA, TRANSB
  INTEGER, INTENT(IN)         :: M, N, K, LDA, LDB, LDC
  REAL(real_8), INTENT(IN)    :: ALPHA, A(LDA,*), B(LDB,*), BETA
  REAL(real_8), INTENT(INOUT) :: C(LDC,*)

  CALL dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
END SUBROUTINE cpmd_dgemm

SUBROUTINE cpmd_dgemmt(UPLO, TRANSA, TRANSB, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  USE kinds,                          ONLY: real_8
  CHARACTER(1), INTENT(IN) :: UPLO, TRANSA, TRANSB
  INTEGER, INTENT(IN) :: N, K, LDA, LDB, LDC
  REAL(real_8), INTENT(IN) :: ALPHA, A(LDA,*), B(LDB,*), BETA
  REAL(real_8), INTENT(INOUT) :: C(LDC,*)
#if defined(_HAS_DGEMMT)
  CALL dgemmt(UPLO, TRANSA, TRANSB, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
#endif
END SUBROUTINE cpmd_dgemmt

SUBROUTINE cpmd_dsymm(SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  USE kinds,                          ONLY: real_8
  CHARACTER(1), INTENT(IN)    :: SIDE, UPLO
  INTEGER, INTENT(IN)         :: M, N, LDA, LDB, LDC
  REAL(real_8), INTENT(IN)    :: ALPHA, A(LDA,*), B(LDB,*), BETA
  REAL(real_8), INTENT(INOUT) :: C(LDC,*)

  CALL dsymm(SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
END SUBROUTINE cpmd_dsymm

SUBROUTINE cpmd_dsyrk(UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC)
  USE kinds,                          ONLY: real_8
  CHARACTER(1), INTENT(IN)    :: UPLO, TRANS
  INTEGER, INTENT(IN)         :: N, K, LDA, LDC
  REAL(real_8), INTENT(IN)    :: ALPHA, A(LDA,*), BETA
  REAL(real_8), INTENT(INOUT) :: C(LDC,*)

  CALL dsyrk(UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC)
END SUBROUTINE cpmd_dsyrk

SUBROUTINE cpmd_dtrmm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
  USE kinds,                          ONLY: real_8
  CHARACTER(1), INTENT(IN)    :: SIDE, UPLO, TRANSA, DIAG
  INTEGER, INTENT(IN)         :: M, N, LDA, LDB
  REAL(real_8), INTENT(IN)    :: ALPHA, A(LDA,*)
  REAL(real_8), INTENT(INOUT) :: B(LDB,*)

  CALL dtrmm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
END SUBROUTINE cpmd_dtrmm

SUBROUTINE cpmd_dger(M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
  USE kinds,                          ONLY: real_8
  INTEGER, INTENT(IN)         :: M, N, INCX, INCY, LDA
  REAL(real_8), INTENT(IN)    :: ALPHA, X(*), Y(*)
  REAL(real_8), INTENT(INOUT) :: A(LDA,*)

  CALL dger(M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
END SUBROUTINE cpmd_dger

SUBROUTINE cpmd_zgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  USE kinds,                          ONLY: real_8
  CHARACTER(1), INTENT(IN)       :: TRANSA, TRANSB
  INTEGER, INTENT(IN)            :: M, N, K, LDA, LDB, LDC
  COMPLEX(real_8), INTENT(IN)    :: ALPHA, A(LDA,*), B(LDB,*), BETA
  COMPLEX(real_8), INTENT(INOUT) :: C(LDC,*)

  CALL zgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
END SUBROUTINE cpmd_zgemm

SUBROUTINE cpmd_zgemv(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
  USE kinds,                          ONLY: real_8
  CHARACTER(1), INTENT(IN)       :: TRANS
  INTEGER, INTENT(IN)            :: M, N, LDA, INCX, INCY
  COMPLEX(real_8), INTENT(IN)    :: ALPHA, X(*), BETA, Y(*)
  COMPLEX(real_8), INTENT(INOUT) :: A(LDA,*)

  CALL zgemv(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
END SUBROUTINE cpmd_zgemv
