      SUBROUTINE EIGVALUE(A, N, E)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Zaikun ZHANG, 26 March, 2013.
C
C     EIGVALUE seeks the eigenvalues of the symmetric matrix A, and 
C     return them in E in ascending order. 
C     Only the upper triangular part of A will be used. 
C     A will not be destroyed. N is the order of A.
C
C     The main computation is done by subroutine DSYEVD of lapack.
C
C     Please contact Zaikun ZHANG (www.zhangzk.net) for any problem
C     concerning EIGVALUE.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         IMPLICIT NONE
         INTEGER (KIND = 4), INTENT(IN) :: N
         REAL (KIND = 8), INTENT(IN) :: A(N, N)
         REAL (KIND = 8), INTENT(OUT) :: E(N)

         REAL (KIND = 8) :: TA(N, N), WORK(1+6*N+2*N*N)
         INTEGER(KIND = 4) :: IWORK(3+5*N), INFO, LDA, LWORK, LIWORK
         CHARACTER :: JOBZ = 'N'
         CHARACTER :: UPLO = 'U'

         LDA = N
         LWORK = 1+6*N+2*N*N
         LIWORK = 3+5*N
         TA = A

         CALL DSYEVD(JOBZ,UPLO, N, TA, LDA, E,WORK, LWORK, IWORK, 
     1    LIWORK, INFO)

         IF (INFO /= 0) WRITE(*,*) "ERROR IN EIGVALUE. DSYEVD 
     1    RETURNS ", INFO

         RETURN
      END
