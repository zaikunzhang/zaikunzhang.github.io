      PROGRAM TEST
        INTEGER (KIND=4) :: N , I, J, K
        REAL (KIND=8) :: D(5)
        REAL (KIND=8) :: B(5,5), G(5), DELTA, L(5), U(5),H(100)
        N = 5
        CALL INIT_RANDOM_SEED()
        CALL RANDOM_NUMBER(B) 
        B = (B+TRANSPOSE(B))/2.0D0
        CALL RANDOM_NUMBER(G) 
        DELTA = 6.0D0
        L = -0.5D0
        U = 0.5D0
        
        K = 0
        DO I = 1, N
            DO J = 1, I
                K = K +1 
                H(K) = B(I,J)
            END DO
        END DO
C        CALL TRS (N, D, B, G, DELTA)
        CALL TRS (N, D, H, G, DELTA)

        WRITE(*,*) D
        WRITE(*,*) SQRT(DOT_PRODUCT(D,D)), DELTA
        WRITE(*,*) G+MATMUL(B,D)
        WRITE(*,*) 0.5D0*DOT_PRODUCT(D,MATMUL(B,D)) + DOT_PRODUCT(D,G)

      END PROGRAM TEST

      SUBROUTINE INIT_RANDOM_SEED()
        INTEGER :: I, N, CLOCK
        INTEGER, DIMENSION(:), ALLOCATABLE :: SEED
      
        CALL RANDOM_SEED(SIZE = N)
        ALLOCATE(SEED(N))
      
        CALL SYSTEM_CLOCK(COUNT=CLOCK)
      
        SEED = CLOCK + 37 * (/ (I - 1, I = 1, N) /)
        CALL RANDOM_SEED(PUT = SEED)
      
        DEALLOCATE(SEED)
      END SUBROUTINE
