/*compare pdgeqpf, pdgeqp3, pdgeqrf, rqrcp
AUTHOR: JIANWEI XIAO AND JULIEN LANGOU*/

/*Reference: 
Fast Parallel Randomized QR with Column Pivoting Algorithms for Reliable Low-rank Matrix Approximations.
Jianwei Xiao, Ming Gu and Julien Langou.
24th IEEE International Conference on High Performance Computing, Data, and Analytics (HIPC), Jaipur, India, 2017.*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

static int max( int a, int b ){
    if (a>b) return(a); else return(b);
}

extern void   Cblacs_pinfo( int* mypnum, int* nprocs);
extern void   Cblacs_get( int context, int request, int* value);
extern int    Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
extern void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
extern void   Cblacs_gridexit( int context);
extern void   Cblacs_exit( int error_code);
extern void Cblacs_barrier(int ConTxt, char *scope);

extern int    numroc_( int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern void   descinit_( int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc,
    int *ictxt, int *lld, int *info);
extern double pdlamch_( int *ictxt , char *cmach);
extern double pdlange_( char *norm, int *m, int *n, double *A, int *ia, int *ja, int *desca, double *work);

extern void pdlacpy_( char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca,
    double *b, int *ib, int *jb, int *descb);
extern void pdgesv_( int *n, int *nrhs, double *A, int *ia, int *ja, int *desca, int* ipiv,
    double *B, int *ib, int *jb, int *descb, int *info);
extern void pdgemm_( char *TRANSA, char *TRANSB, int * M, int * N, int * K, double * ALPHA,
    double * A, int * IA, int * JA, int * DESCA, double * B, int * IB, int * JB, int * DESCB,
    double * BETA, double * C, int * IC, int * JC, int * DESCC );
extern int  indxg2p_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);

extern int rqrcp_( int *m, int *n, int *k, double *A, int *descA, int *m_B, int *n_B, double *B, int *descB, double *OMEGA, int *desc_OMEGA, int *ipiv, double *tau, int *nb,
    int *ipiv_a, double *tau_b, double *work, int *lwork );

extern double check_qpf_( int *m, int *n, int *k, double *A, int *descA, double *B, int *descB, double *C, int *descC, 
    double *tau, int *ipiv, double *work, int *lwork );

extern double check_qrf_( int *m, int *n, int *k, double *A, int *descA, double *B, int *descB, double *C, int *descC, 
    double *tau, double *work, int *lwork );

extern void  partial_pdgeqpf_( int *k, int *m, int *n, double *a, int *ia, int *ja, int *desca, int *ipiv, double *tau,
    double *work, int *lwork, int *info );

extern void partial_pdgeqrf_ ( int *k, int *m, int *n, double *a, int *ia, int *ja, int *descA, double *tau, double *work, int *lwork, int *info);

extern void  partial_pdgeqp3_( int *k, int *m, int *n, double *a, int *desca, int *ipiv, double *tau,
    double *work, int *lwork, int *iwork, int *liwork, int *info );

double randn(double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}

int main(int argc, char **argv) {
    int iam, nprocs;
    int myrank_mpi, nprocs_mpi;
    int ictxt, nprow, npcol, myrow, mycol;
    int mp, nq, m, n, k, nb_dist, nb_alg, np_oversampling;
    int ii, i, j, info, itemp, seed, lwork, lwork_save, liwork;
    int descA[9];
    double *A, *Acpy, *work, *tau, *C, *Origin;
    int *ipiv, *iwork;
    int izero=0,ione=1;

    double bwderr;

    double MPIt1, MPIt2, MPIelapsed, GFLOPS, duration, global;

    MPI_Init( &argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);

    int len;
    char name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(name, &len);
    //printf("Hello, world.  I am %d of %d on %s\n", myrank_mpi, nprocs_mpi, name);
    fflush(stdout); 

    m = 500; n = 2000; k = 500; nprow = 2; npcol = 3; nb_alg = 32; np_oversampling = 10; nb_dist = 256;

    for( i = 1; i < argc; i++ ) {
        if( strcmp( argv[i], "-m" ) == 0 ) {
            m = atoi(argv[i+1]);
            i++;
        }
        if( strcmp( argv[i], "-n" ) == 0 ) {
            n = atoi(argv[i+1]);
            i++;
        }
        if( strcmp( argv[i], "-k" ) == 0 ) {
            k = atoi(argv[i+1]);
            i++;
        }
        if( strcmp( argv[i], "-p" ) == 0 ) {
            nprow = atoi(argv[i+1]);
            i++;
        }
        if( strcmp( argv[i], "-q" ) == 0 ) {
            npcol = atoi(argv[i+1]);
            i++;
        }
        if( strcmp( argv[i], "-nb_alg" ) == 0 ) {
            nb_alg = atoi(argv[i+1]);
            i++;
        }
        if( strcmp( argv[i], "-np_oversampling" ) == 0 ) {
            np_oversampling = atoi(argv[i+1]);
            i++;
        }
        if( strcmp( argv[i], "-nb_dist" ) == 0 ) {
            nb_dist = atoi(argv[i+1]);
            i++;
        }
    }

    if (nprow*npcol>nprocs_mpi){
        if (myrank_mpi==0){
            printf(" **** ERROR : we do not have enough processes available to make a p-by-q process grid ***\n");
            printf(" **** Bye-bye                                                                         ***\n");
        }
        MPI_Finalize(); exit(1);
    }

    Cblacs_pinfo( &iam, &nprocs ) ;
    Cblacs_get( -1, 0, &ictxt );
    Cblacs_gridinit( &ictxt, "Row", nprow, npcol );
    Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

    if ((myrow < nprow)&(mycol < npcol)){

        mp = numroc_( &m, &nb_dist, &myrow, &izero, &nprow );
        nq = numroc_( &n, &nb_dist, &mycol, &izero, &npcol );

        seed = iam*m*n; srand(seed);

        A = (double *)malloc(mp*nq*sizeof(double)) ;
        Acpy = (double *)malloc(mp*nq*sizeof(double)) ;
        Origin = (double *)malloc(mp*nq*sizeof(double)) ;

        ii = 0;
        for (i = 0; i < mp; i++) {
            for (j = 0; j < nq; j++) {
                Origin[ii] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
                ii++;   
            }
        }

        itemp = max( 1, mp );
        descinit_( descA, &m, &n, &nb_dist, &nb_dist, &izero, &izero, &ictxt, &itemp, &info );

/*      PDGEQRF */
        pdlacpy_( "All", &m, &n, Origin, &ione, &ione, descA, A, &ione, &ione, descA );
        pdlacpy_( "All", &m, &n, A, &ione, &ione, descA, Acpy, &ione, &ione, descA );

        tau = (double *)malloc((nq+nb_dist)*sizeof(double)) ;

        lwork = -1;
        work = (double *)malloc(1*sizeof(double)) ;
        partial_pdgeqrf_( &k, &m, &n, A, &ione, &ione, descA, tau, work, &lwork, &info );
        lwork = (int) work[0];
        free(work);
        work = (double *)malloc(lwork*sizeof(double)) ;

        lwork_save = lwork;

        pdlacpy_( "All", &m, &n, Origin, &ione, &ione, descA, A, &ione, &ione, descA );
        pdlacpy_( "All", &m, &n, A, &ione, &ione, descA, Acpy, &ione, &ione, descA );

        Cblacs_barrier( ictxt, "A" );
        MPIt1 = MPI_Wtime();
        partial_pdgeqrf_( &k, &m, &n, A, &ione, &ione, descA, tau, work, &lwork, &info );
        Cblacs_barrier( ictxt, "A" );
        MPIt2 = MPI_Wtime();
        MPIelapsed=MPIt2-MPIt1;
        MPI_Reduce(&MPIelapsed,&global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
        GFLOPS=4.0e+00/3.0e+00*((double) m)*((double) n)*((double) k)/MPIelapsed*1.0e-9;

        C = (double *)malloc(mp*nq*sizeof(double)) ;
        bwderr = check_qrf_( &m, &n, &k, Acpy, descA, A, descA, C, descA, tau, work, &lwork );
        free( C );

        if ( iam==0 ){
            printf("\t\t nprow \t npcol \t m \t n \t k \t nb_dist nb_alg \t\t correctness(zero) MPIelapsed GFLOPS");
            printf("\n");
        }

        if ( iam==0 ){
            printf("PDGEQRF \t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t     \t     \t%e\t%f\t%f\n", nprow, npcol, m, n, k, nb_dist, nb_alg, bwderr, global, GFLOPS);
        }

        free(work);
        free(tau);

/*      rqrcp     */
        pdlacpy_( "All", &m, &n, Origin, &ione, &ione, descA, A, &ione, &ione, descA );
        pdlacpy_( "All", &m, &n, A, &ione, &ione, descA, Acpy, &ione, &ione, descA );

        {
            int m_omega = nb_alg + np_oversampling;
            int mp_omega, mq;
            double *OMEGA, *B, *tau_B;
            int descB[9], desc_OMEGA[9];
            int *ipiv_B;

            mp_omega = numroc_( &m_omega, &nb_dist, &myrow, &izero, &nprow );
            mq = numroc_( &m, &nb_dist, &mycol, &izero, &npcol );

            OMEGA = (double *)malloc(mp_omega*mq*sizeof(double)) ;
            B = (double *)malloc(mp_omega*nq*sizeof(double)) ;

            tau = (double *)malloc((nq+nb_dist)*sizeof(double)) ;
            ipiv = (int *)malloc((nq+nb_dist)*sizeof(int)) ;
            tau_B = (double *)malloc((nq+nb_dist)*sizeof(double)) ;
            ipiv_B = (int *)malloc((nq+nb_dist)*sizeof(int)) ;

            lwork = lwork_save;
            work = (double *)malloc(lwork*sizeof(double)) ;

            ii = 0;
            for (i = 0; i < mp_omega; i++) {
                for (j = 0; j < mq; j++) {
                    OMEGA[ii] = randn(0.0, 1.0);
                    ii++;   
                }
            }

            itemp = max( 1, mp_omega ); descinit_( descB, &m_omega, &n, &nb_dist, &nb_dist, &izero, &izero, &ictxt, &itemp, &info );
            itemp = max( 1, mp_omega ); descinit_( desc_OMEGA, &m_omega, &m, &nb_dist, &nb_dist, &izero, &izero, &ictxt, &itemp, &info );

            Cblacs_barrier( ictxt, "A" );
            MPIt1 = MPI_Wtime();

            rqrcp_( &m, &n, &k, A, descA, &m_omega, &n, B, descB, OMEGA, desc_OMEGA, ipiv, tau, &nb_alg, ipiv_B, tau_B, work, &lwork );

            Cblacs_barrier( ictxt, "A" );
            MPIt2 = MPI_Wtime();
            MPIelapsed=MPIt2-MPIt1;
            MPI_Reduce(&MPIelapsed,&global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
            GFLOPS=4.0e+00/3.0e+00*((double) m)*((double) n)*((double) k)/MPIelapsed*1.0e-9;

            C = (double *)malloc(mp*nq*sizeof(double)) ;
            bwderr = check_qpf_( &m, &n, &k, Acpy, descA, A, descA, C, descA, tau, ipiv, work, &lwork );
            free(work);
            free( C );

            if ( iam==0 ){
                printf("RQRCP  \t\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t     \t     \t%e\t%f\t%f\n", nprow, npcol, m, n, k, nb_dist, nb_alg, bwderr, global, GFLOPS);
            }

            free( OMEGA );
            free( B );
            free(tau);
            free(ipiv);
            free(tau_B);
            free(ipiv_B);

        }

/*      PDGEQP3 */

        pdlacpy_( "All", &m, &n, Origin, &ione, &ione, descA, A, &ione, &ione, descA );
        pdlacpy_( "All", &m, &n, A, &ione, &ione, descA, Acpy, &ione, &ione, descA );

        tau = (double *)malloc((nq+nb_dist)*sizeof(double)) ;
        ipiv = (int *)malloc((nq+nb_dist)*sizeof(int)) ;

        lwork = 2*lwork_save;
        liwork = 2*lwork_save;
        work = (double *)malloc(lwork*sizeof(double)) ;
        iwork = (int *)malloc(liwork*sizeof(int));

        pdlacpy_( "All", &m, &n, Origin, &ione, &ione, descA, A, &ione, &ione, descA );
        pdlacpy_( "All", &m, &n, A, &ione, &ione, descA, Acpy, &ione, &ione, descA );

        Cblacs_barrier( ictxt, "A" );
        MPIt1 = MPI_Wtime();
        partial_pdgeqp3_( &k, &m, &n, A, descA, ipiv, tau, work, &lwork, iwork, &liwork, &info );
        Cblacs_barrier( ictxt, "A" );
        MPIt2 = MPI_Wtime();
        MPIelapsed=MPIt2-MPIt1;
        MPI_Reduce(&MPIelapsed,&global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
        GFLOPS=4.0e+00/3.0e+00*((double) m)*((double) n)*((double) k)/MPIelapsed*1.0e-9;

        C = (double *)malloc(mp*nq*sizeof(double)) ;
        lwork = lwork_save;
        free(work);
        work = (double *)malloc(lwork*sizeof(double)) ;
        bwderr = check_qpf_( &m, &n, &k, Acpy, descA, A, descA, C, descA, tau, ipiv, work, &lwork );
        free(work);
        free( C );

        if ( iam==0 ){
            printf("PDGEQP3 \t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t     \t     \t%e\t%f\t%f\n", nprow, npcol, m, n, k, nb_dist, nb_alg, bwderr, global, GFLOPS);
        }

        free(iwork);
        free(tau);
        free(ipiv);


/*      PDGEQPF */

        pdlacpy_( "All", &m, &n, Origin, &ione, &ione, descA, A, &ione, &ione, descA );
        pdlacpy_( "All", &m, &n, A, &ione, &ione, descA, Acpy, &ione, &ione, descA );

        tau = (double *)malloc((nq+nb_dist)*sizeof(double)) ;
        ipiv = (int *)malloc((nq+nb_dist)*sizeof(int)) ;

        lwork = 2 * lwork_save;
        work = (double *)malloc(lwork*sizeof(double)) ;

        pdlacpy_( "All", &m, &n, Origin, &ione, &ione, descA, A, &ione, &ione, descA );
        pdlacpy_( "All", &m, &n, A, &ione, &ione, descA, Acpy, &ione, &ione, descA );

        Cblacs_barrier( ictxt, "A" );
        MPIt1 = MPI_Wtime();
        partial_pdgeqpf_( &k, &m, &n, A, &ione, &ione, descA, ipiv, tau, work, &lwork, &info );
        Cblacs_barrier( ictxt, "A" );
        MPIt2 = MPI_Wtime();
        MPIelapsed=MPIt2-MPIt1;
        MPI_Reduce(&MPIelapsed,&global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
        GFLOPS=4.0e+00/3.0e+00*((double) m)*((double) n)*((double) k)/MPIelapsed*1.0e-9;

        C = (double *)malloc(mp*nq*sizeof(double)) ;
        lwork = lwork_save;
        free(work);
        work = (double *)malloc(lwork*sizeof(double)) ;
        bwderr = check_qpf_( &m, &n, &k, Acpy, descA, A, descA, C, descA, tau, ipiv, work, &lwork );
        free( C );

        if ( iam==0 ){
            printf("PDGEQPF \t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t     \t     \t%e\t%f\t%f\n", nprow, npcol, m, n, k, nb_dist, nb_alg, bwderr, global, GFLOPS);
        }

        free(work);
        free(tau);
        free(ipiv);

/*                 */
        free(A);
        free(Acpy);

    }

    Cblacs_gridexit( 0 );
    MPI_Finalize();
    exit(0);
}

