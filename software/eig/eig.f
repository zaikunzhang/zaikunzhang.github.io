      SUBROUTINE EIG(A, N, K, LAMBDA, V)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Zaikun ZHANG, 26 March, 2013.
C
C     EIG seeks the largest K eigenvalues of A, and return them in
C     LAMBDA, in descending order, and the columns of V are 
C     corresponding eigenvectors. A is a symmetric matrix of order N. 
C     Only the upper triangular part of A will be used. 
C     A will not be destroyed.
C
C     The main computation is done by subroutine DSYEVD of lapack.
C
C     Please contact Zaikun ZHANG (www.zhangzk.net) for any problem
C     concerning EIG.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         IMPLICIT NONE

         INTEGER (KIND = 4), INTENT(IN) :: N, K
         REAL(KIND = 8), INTENT(IN) :: A(N, N) 
         REAL(KIND = 8), INTENT(OUT) :: LAMBDA(K), V(N, K) 

         REAL(KIND = 8) :: WORK(1+6*N+2*N**2), TEMPA(N, N), W(N)
         INTEGER(KIND = 4) :: IWORK(3+5*N) 
         INTEGER(KIND = 4) :: LWORK, LIWORK, INFO, I
         CHARACTER :: JOBZ, UPLO

         JOBZ = 'V'
         UPLO = 'U'
         LWORK = 1+6*N+2*N**2
         LIWORK = 3+5*N
         TEMPA = A

         CALL DSYEVD(JOBZ, UPLO, N, TEMPA, N, W, WORK, LWORK, IWORK,
     1   LIWORK, INFO)

         IF (INFO /= 0) WRITE(*,*) "ERROR IN EIG. DSYEVD RETURNS ", INFO
         
         DO I = 1, K
            LAMBDA(I) = W(N-I+1)
            V(:,I) = TEMPA(:, N-I+1)
         END DO

         RETURN
      END
