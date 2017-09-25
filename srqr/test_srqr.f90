! Reference: 
! Fast Parallel Randomized QR with Column Pivoting Algorithms for Reliable Low-rank Matrix Approximations.
! Jianwei Xiao, Ming Gu and Julien Langou.
! 24th IEEE International Conference on High Performance Computing, Data, and Analytics (HIPC), Jaipur, India, 2017.

PROGRAM TEST_RQRCP_SWAP

	IMPLICIT NONE 

	CHARACTER(10):: M_STRING, N_STRING, K_STRING, NB_MAX_STRING, P_STRING, G_STRING, D_STRING, COUNTER_STRING
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A, A_SRQR, A_DGEQPF, A_DGEQRF, HOUSEHOLDER, GIVENS_ARRAY, EXTRA_HOUSEHOLDER
	INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TAU, WORK, EXTRA_TAU
	DOUBLE PRECISION :: DLANGE, A_NORM, C, S, SCALE, DTEMP, G, START, END
	INTEGER :: ISEED(4) = (/38,56,67,86/), I, J, INFO, COUNTER, LWORK, M, N, K, NB_MAX, P, D, NUM_SWAPS, UPPER_BOUND
	INTEGER, DIMENSION(:), ALLOCATABLE :: INDEX_ARRAY
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AK
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AKS1, AKS2, AKS3, AS
	INTEGER :: PRINT_CHOICE
	DOUBLE PRECISION :: U, VT

	CALL GETARG(1, M_STRING)
	CALL GETARG(2, N_STRING)
	CALL GETARG(3, K_STRING)
	CALL GETARG(4, NB_MAX_STRING)
	CALL GETARG(5, P_STRING)
	CALL GETARG(6, G_STRING)
	CALL GETARG(7, D_STRING)
	CALL GETARG(8, COUNTER_STRING)

	READ(M_STRING, *) M
	READ(N_STRING, *) N
	READ(K_STRING, *) K
	READ(NB_MAX_STRING, *) NB_MAX
	READ(P_STRING, *) P
	READ(G_STRING, *) G
	READ(D_STRING, *) D
	READ(COUNTER_STRING, *) COUNTER

	LWORK = MAX(1,128*N,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
	UPPER_BOUND = 20
	PRINT_CHOICE = 1
	ALLOCATE(AKS1(K), AKS2(K), AKS3(K), AS(MIN(M,N)))
	ALLOCATE(A(M,N), A_SRQR(M,N), A_DGEQPF(M,N), A_DGEQRF(M,N), IPIV(N), TAU(N), WORK(LWORK), AK(K,K))
	ALLOCATE(HOUSEHOLDER(M,K+1), GIVENS_ARRAY(K,UPPER_BOUND), INDEX_ARRAY(UPPER_BOUND), EXTRA_HOUSEHOLDER(M-K,UPPER_BOUND), &
		EXTRA_TAU(UPPER_BOUND))

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF (COUNTER .EQ. 0) THEN
		A = 0.0D+0
		A(1,1) = 1.0D+0
		! CALL DLARNV(3, ISEED, M * N, A)
		IF (PRINT_CHOICE .EQ. 1) THEN
			WRITE (*,*) 'TEST ON RANDOM MATRIX.'
			WRITE (*,*) 'M =', M
			WRITE (*,*) 'N =', N
			WRITE (*,*) 'K =', K
			WRITE (*,*) 'NB_MAX =', NB_MAX
			WRITE (*,*) 'P =', P
			WRITE (*,*) 'G =', G
			WRITE (*,*) 'D =', D
			WRITE (*,*)
		ELSE
		END IF
	ELSE
	END IF

	IF (COUNTER .EQ. 1) THEN
		OPEN(UNIT = 1, FILE = "HAPT_7767_561.DAT", ACTION = "READ")
		DO J = 1,N
			DO I = 1,M
				READ(1,*) A(I,J)
			END DO
		END DO
		IF (PRINT_CHOICE .EQ. 1) THEN
			WRITE (*,*) 'TEST ON HAPT, 7767 * 561.'
			WRITE (*,*) 'M =', M
			WRITE (*,*) 'N =', N
			WRITE (*,*) 'K =', K
			WRITE (*,*) 'NB_MAX =', NB_MAX
			WRITE (*,*) 'P =', P
			WRITE (*,*) 'G =', G
			WRITE (*,*) 'D =', D
			WRITE (*,*)
		ELSE
		END IF
	ELSE
	END IF

	IF (COUNTER .EQ. 2) THEN
		OPEN(UNIT = 1, FILE = "MNIST_60000_784.DAT", ACTION = "READ")
		DO J = 1,N
			DO I = 1,M
				READ(1,*) A(I,J)
			END DO
		END DO
		IF (PRINT_CHOICE .EQ. 1) THEN
			WRITE (*,*) 'TEST ON MNIST, 60000 * 784.'
			WRITE (*,*) 'M =', M
			WRITE (*,*) 'N =', N
			WRITE (*,*) 'K =', K
			WRITE (*,*) 'NB_MAX =', NB_MAX
			WRITE (*,*) 'P =', P
			WRITE (*,*) 'G =', G
			WRITE (*,*) 'D =', D
			WRITE (*,*)
		ELSE
		END IF
	ELSE
	END IF

	IF (COUNTER .EQ. 3) THEN
		C = 0.285D+0
		! S = SQRT((1.0D+0-C) * (1.0D+0+C))
		S = SQRT((0.9999D+0-C) * (0.9999D+0+C))
		SCALE = 1.0D+0
		A = 0.0D+0
		DO I = 1,N
			A(I,I) = 1.0D+0
		END DO

		DO I = 1,N
			DO J = I+1,N
				A(I,J) = -C
			END DO
		END DO

		DO I = 1,N
			DO J = 1,N
				A(I,J) = A(I,J) * SCALE
			END DO
			SCALE = SCALE * S
		END DO

		IF (PRINT_CHOICE .EQ. 1) THEN
			WRITE (*,*) 'TEST ON A KAHAN MATRIX'
			WRITE (*,*) 'M =', M
			WRITE (*,*) 'N =', N
			WRITE (*,*) 'K =', K
			WRITE (*,*) 'NB_MAX =', NB_MAX
			WRITE (*,*) 'P =', P
			WRITE (*,*) 'G =', G
			WRITE (*,*) 'D =', D
			WRITE (*,*)
		ELSE
		END IF
	ELSE
	END IF

	CALL DLACPY('ALL', M, N, A, M, A_SRQR, M)
	CALL DLACPY('ALL', M, N, A, M, A_DGEQPF, M)
	CALL DLACPY('ALL', M, N, A, M, A_DGEQRF, M)
	A_NORM = DLANGE('F', M, N, A, M, WORK)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! RUN SRQR ON A
	! WRITE (*,*)
	IF (PRINT_CHOICE .EQ. 1) THEN
		WRITE (*,*) 'SRQR'
	ELSE
	END IF
	IPIV = 0
	TAU = 0.0D+0
	CALL CPU_TIME(START)
	CALL SRQR(M, N, K, NB_MAX, P, A_SRQR, TAU, IPIV, G, D, WORK, LWORK, NUM_SWAPS, UPPER_BOUND, HOUSEHOLDER, &
		GIVENS_ARRAY, INDEX_ARRAY, EXTRA_HOUSEHOLDER, EXTRA_TAU)
	CALL CPU_TIME(END)
	IF (PRINT_CHOICE .EQ. 1) THEN
		WRITE (*,*) 'RUN TIME IS', END - START
		WRITE (*,*) 'NUMBER OF EXTRA SWAPS', NUM_SWAPS
		WRITE (*,*) 'RESIDUAL_NORM / A_NORM, MEASURE APPROXIMATION QUALITY: ', DLANGE('F', M-K, N-K, A_SRQR(K+1,K+1), M, WORK) / A_NORM
		WRITE (*,*) 
	ELSE
		WRITE (*,*) END - START
		WRITE (*,*) DLANGE('F', M-K, N-K, A_SRQR(K+1,K+1), M, WORK) / A_NORM
	END IF
	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! RUN DGEQPF ON A_DGEQPF
	! WRITE (*,*)
	IF (PRINT_CHOICE .EQ. 1) THEN
		WRITE (*,*) 'QRCP'
	ELSE
	END IF
	IPIV = 0
	TAU = 0.0D+0
	CALL CPU_TIME(START)
	CALL PARTIAL_DGEQPF(K, M, N, A_DGEQPF, M, IPIV, TAU, WORK, INFO)
	CALL CPU_TIME(END)
	IF (PRINT_CHOICE .EQ. 1) THEN
		WRITE (*,*) 'RUN TIME IS', END - START
		WRITE (*,*) 'RESIDUAL_NORM / A_NORM, MEASURE APPROXIMATION QUALITY: ', DLANGE('F', M-K, N-K, A_DGEQPF(K+1,K+1), M, WORK) / A_NORM
		WRITE (*,*) 
	ELSE
		WRITE (*,*) END - START
		WRITE (*,*) DLANGE('F', M-K, N-K, A_DGEQPF(K+1,K+1), M, WORK) / A_NORM
	END IF
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! RUN DGEQRF ON A_DGEQRF
	! WRITE (*,*)
	IF (PRINT_CHOICE .EQ. 1) THEN
		WRITE (*,*) 'UNPIVOTED QR'
	ELSE
	END IF

	IF (COUNTER .NE. 3) THEN
		DO J = 1,N
			IPIV(J) = J
		END DO
		CALL RANDPERM(N, IPIV)
		CALL DLAPMT('FORWARD', M, N, A_DGEQRF, M, IPIV)
	ELSE
	END IF
	IPIV = 0
	TAU = 0.0D+0

	CALL CPU_TIME(START)
	CALL DGEQRF(M, K, A_DGEQRF, M, TAU, WORK, LWORK, INFO)
	CALL DORMQR('L', 'T', M, N-K, K, A_DGEQRF, M, TAU, A_DGEQRF(1,K+1), M, WORK, LWORK, INFO) 
	CALL CPU_TIME(END)
	IF (PRINT_CHOICE .EQ. 1) THEN
		WRITE (*,*) 'RUN TIME IS', END - START
		WRITE (*,*) 'RESIDUAL_NORM / A_NORM, MEASURE APPROXIMATION QUALITY: ', DLANGE('F', M-K, N-K, A_DGEQRF(K+1,K+1), M, WORK) / A_NORM
		WRITE (*,*) 
	ELSE
		WRITE (*,*) END - START
		WRITE (*,*) DLANGE('F', M-K, N-K, A_DGEQRF(K+1,K+1), M, WORK) / A_NORM
	END IF

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! PRINT THE SINGULAR VALUE RATIO 
	! ONLY PRINT FOR KAHAN MATRIX
	IF (COUNTER .EQ. 3) THEN
		CALL DGESVD('N', 'N', M, N, A, M, AS, U, M, VT, N, WORK, LWORK, INFO)

		CALL DLASET('ALL', K, K, 0.0D+0, 0.0D+0, AK, K)
		CALL DLACPY('U', K, K, A_SRQR, M, AK, K)
		CALL DGESVD('N', 'N', K, K, AK, K, AKS1, U, K, VT, K, WORK, LWORK, INFO)

		CALL DLASET('ALL', K, K, 0.0D+0, 0.0D+0, AK, K)
		CALL DLACPY('U', K, K, A_DGEQPF, M, AK, K)
		CALL DGESVD('N', 'N', K, K, AK, K, AKS2, U, K, VT, K, WORK, LWORK, INFO)

		CALL DLASET('ALL', K, K, 0.0D+0, 0.0D+0, AK, K)
		CALL DLACPY('U', K, K, A_DGEQRF, M, AK, K)
		CALL DGESVD('N', 'N', K, K, AK, K, AKS3, U, K, VT, K, WORK, LWORK, INFO)

		WRITE (*,*) 'SIGMA_J(R_11)/SIGMA_J(A)'
		WRITE (*,'(20G12.4)') '     INDEX ', 'SRQR ', 'QRCP ', 'QR'
		DO I = 1,K
		   WRITE (*,'(20G12.4)') I, AKS1(I)/AS(I), AKS2(I)/AS(I), AKS3(I)/AS(I)
		END DO
		WRITE (*,*)
	ELSE
	END IF

	DEALLOCATE(A, A_SRQR, A_DGEQPF, A_DGEQRF, IPIV, TAU, WORK, AK)
	DEALLOCATE(HOUSEHOLDER, GIVENS_ARRAY, INDEX_ARRAY)
	DEALLOCATE(AKS1, AKS2, AKS3, AS)
END PROGRAM TEST_RQRCP_SWAP
