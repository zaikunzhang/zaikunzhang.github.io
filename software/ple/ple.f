      SUBROUTINE PLE(N, A, b)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Zaikun ZHANG, 26 March, 2013.
C
C     PLE solves the linear equation Ax=b, where A is a positive
C     definite matrix of order N. Only the upper triangular part of A
C     will be used.
C     A will be destroyed. The solution will be returned in b.
C
C     The main computation is done by subroutine DPOTRF and DPOTRS
C     of lapack.
C
C     Please contact Zaikun ZHANG (www.zhangzk.net) for any problem
C     concerning PLE.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        INTEGER(KIND=4), INTENT(IN) :: N
        REAL(KIND=8), INTENT(INOUT) :: A(N,N), b(N)

        CHARACTER UPLO 
        INTEGER(KIND=4) :: LDA, LDB, NRHS
        INTEGER(KIND=4) :: INFO

        UPLO = 'U'
        LDA = N
        LDB = N
        NRHS = 1
        
        CALL DPOTRF(UPLO, N, A, LDA, INFO)
        IF(INFO /= 0) THEN
            WRITE(*,*) "ERROR IN PLE. DPOTRF EXITS WITH INFO = ", INFO
        END IF

        CALL DPOTRS(UPLO, N, NRHS, A, LDA, b, LDB, INFO)
        IF(INFO /= 0) THEN
            WRITE(*,*) "ERROR IN PLE. DPOTRS EXITS WITH INFO = ", INFO
        END IF

        RETURN
      END SUBROUTINE PLE
