      SUBROUTINE TRB (N, D, B, G, DELTA, L, U)
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Zaikun ZHANG, 26 March, 2013.
C
C     TRB seeks an inexact solution to the trust region subproblem 
C     with box constraints: 
C     min 0.5*D'*B*D + G'*D
C     s.t. ||D|| <= DELTA,
C          L <= D <= U.
C     Only the lower triangular part of B is used. 
C     The main computation is done by subroutine TRSBOX in BOBYQA 
C     (bobyqa.f) by Professor M.J.D. Powell. Please see line 68 of 
C     this subroutine and lines 1668 -- 2068 of bobyqa.f for details. 
C
C     TRB is only applicable when the trust region center 0 is feasible.
C
C     Please ask Professor Powell (mjdp@damtp.cam.ac.uk) for the 
C     code of bobyqa.f, and see 
C     www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf
C     for a paper on BOBYQA.
C
C     Please contact Zaikun ZHANG (www.zhangzk.net) for any problem
C     concerning TRB.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        IMPLICIT NONE

CCCCCCCCCCCCCCCCCCCC BEGIN DUMMIES CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        INTEGER (KIND=4), INTENT (IN) :: N
        REAL (KIND=8), INTENT (OUT) :: D(N)
        REAL (KIND=8), INTENT (IN) :: B(N,N), G(N), DELTA, L(N), U(N)
CCCCCCCCCCCCCCCCCCCC END DUMMIES CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCC BEGIN WORKING-VARIABLES CCCCCCCCCCCCCCCCCCCCCCCCCC
        INTEGER (KIND=4) :: NPT, I, J, IH
        REAL (KIND=8) :: XPT(0,N), XOPT(N), HQ(N*(N+1)/2), PQ(0), 
     1   XNEW(N), GNEW(N), XBDI(N), S(N), HS(N), HRED(N), DSQ, CRVMIN
CCCCCCCCCCCCCCCCCCCC END WORKING-VARIABLES CCCCCCCCCCCCCCCCCCCCCCCCCCCC
        
CCCCCCCCCCCCCCCCCCCC BEGIN MAIN-PROCEDURE CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        IF (MAXVAL(L) > 0.0D0 .OR. MINVAL(U) < 0.0D0) THEN
            WRITE(*,*) "TRB returns because the trust region center is
     1       infeasible."
            RETURN
        END IF

        NPT = 0
        XPT = 0.0D0
        XOPT = 0.0D0
        IH = 0
        DO I = 1, N
          DO J = 1, I
              IH = IH + 1
              HQ(IH) = B(I, J)
          END DO
        END DO
        PQ = 0.0D0
        XNEW = 0.0D0
        D = 0.0D0
        GNEW = 0.0D0
        XBDI = 0.0D0
        S = 0.0D0
        HRED = 0.0D0
        DSQ = 0.0D0
        CRVMIN = 0.0D0

        CALL TRSBOX (N, NPT, XPT, XOPT, G, HQ, PQ, L, U, DELTA,
     1   XNEW, D, GNEW, XBDI, S, HS, HRED, DSQ, CRVMIN)
CCCCCCCCCCCCCCCCCCCC END MAIN-PROCEDURE CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        RETURN

      END SUBROUTINE TRB


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC     Another version. The only difference lies in the data structure of 
CC     the matrix B.
C      SUBROUTINE TRB (N, D, B, G, DELTA, L, U)
C      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC     Zaikun ZHANG, 26 March, 2013.
CC
CC     TRB seeks an inexact solution to the trust region subproblem 
CC     with box constraints: 
CC     min 0.5*D'*B*D + G'*D
CC     s.t. ||D|| <= DELTA,
CC          L <= D <= U.
CC     B holds the lower triangular part of B, in the order B(1,1),  
CC     B(2,1), B(2,2), B(3,1), B(3,2), B(3,3), ...   
CC     The main computation is done by subroutine TRSBOX in BOBYQA 
CC     (bobyqa.f) by Professor M.J.D. Powell. Please see line 62 of 
CC     this subroutine and lines 1668 -- 2068 of bobyqa.f for details. 
CC
CC     TRB is only applicable when the trust region center 0 is feasible.
CC
CC     Please ask Professor Powell (mjdp@damtp.cam.ac.uk) for the 
CC     code of bobyqa.f, and see 
CC     www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf
CC     for a paper on BOBYQA.
CC
CC     Please contact Zaikun ZHANG (www.zhangzk.net) for any problem
CC     concerning TRB.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        IMPLICIT NONE
C
CCCCCCCCCCCCCCCCCCCCC BEGIN DUMMIES CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        INTEGER (KIND=4), INTENT (IN) :: N
C        REAL (KIND=8), INTENT (OUT) :: D(N)
C        REAL (KIND=8), INTENT (IN) :: B(N,N), G(N), DELTA, L(N), U(N)
CCCCCCCCCCCCCCCCCCCCC END DUMMIES CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCC BEGIN WORKING-VARIABLES CCCCCCCCCCCCCCCCCCCCCCCCCC
C        INTEGER (KIND=4) :: NPT
C        REAL (KIND=8) :: XPT(0,N), XOPT(N), PQ(0), 
C     1   XNEW(N), GNEW(N), XBDI(N), S(N), HS(N), HRED(N), DSQ, CRVMIN
CCCCCCCCCCCCCCCCCCCCC END WORKING-VARIABLES CCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        
CCCCCCCCCCCCCCCCCCCCC BEGIN MAIN-PROCEDURE CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        IF (MAXVAL(L) > 0.0D0 .OR. MINVAL(U) < 0.0D0) THEN
C            WRITE(*,*) "TRB returns because the trust region center is
C     1       infeasible."
C            RETURN
C        END IF
C
C        NPT = 0
C        XPT = 0.0D0
C        XOPT = 0.0D0
C        PQ = 0.0D0
C        XNEW = 0.0D0
C        D = 0.0D0
C        GNEW = 0.0D0
C        XBDI = 0.0D0
C        S = 0.0D0
C        HRED = 0.0D0
C        DSQ = 0.0D0
C        CRVMIN = 0.0D0
C
C        CALL TRSBOX (N, NPT, XPT, XOPT, G, B, PQ, L, U, DELTA,
C     1   XNEW, D, GNEW, XBDI, S, HS, HRED, DSQ, CRVMIN)
CCCCCCCCCCCCCCCCCCCCC END MAIN-PROCEDURE CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        RETURN
C
C      END SUBROUTINE TRB
