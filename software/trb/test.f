      PROGRAM TEST
        INTEGER (KIND=4) :: N 
        REAL (KIND=8) :: D(5)
        REAL (KIND=8) :: B(5,5), G(5), DELTA, L(5), U(5)
        N = 5
        CALL INIT_RANDOM_SEED()
        CALL RANDOM_NUMBER(B) 
        B = (B+TRANSPOSE(B))/2.0D0
        CALL RANDOM_NUMBER(G) 
        DELTA = 0.6D0
        L = -0.5D0
        U = 0.5D0

        CALL TRB (N, D, B, G, DELTA, L, U)

        WRITE(*,*) D
        WRITE(*,*) SQRT(DOT_PRODUCT(D,D)), DELTA
        WRITE(*,*) D-L
        WRITE(*,*) U-D
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
