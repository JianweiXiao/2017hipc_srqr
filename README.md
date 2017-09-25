# 2017hipc_srqr

These folders contain the codes for numerical experiments in the paper 
Fast Parallel Randomized QR with Column Pivoting Algorithms for Reliable Low-rank Matrix Approximations.
Jianwei Xiao, Ming Gu and Julien Langou.
24th IEEE International Conference on High Performance Computing, Data, and Analytics (HIPC), Jaipur, India, 2017.

prerequisites: LAPACK, ScaLAPACK

notice: you may need to modify make.inc to be suitable for your system

*********************************************************************************************

Folder “srqr” contains the fortran code for numerical experiments A and B. (sequential srqr)

Implementation:
enter make, then ./job_compare_hapt.sh ./job_compare_kahan.sh
mnist data required by ./job_compare_mnist.sh is not included because it's too large to upload to github

hapt data comes from
J.-L. Reyes-Ortiz, L. Oneto, A. Sama, X. Parra, and D. An-guita,  “Transition-aware  human  activity  recognition  usingsmartphones,”Neurocomputing, vol. 171, pp. 754–767, 2016.

*********************************************************************************************

Folder “dist_rqrcp” contains the fortran code for numerical experiment C. (parallel rqrcp)

Implementation:
enter make, then ./job_xtimer_qr.sh

*********************************************************************************************

Citing this work:

We ask those who benefit from this work to cite the paper 
Fast Parallel Randomized QR with Column Pivoting Algorithms for Reliable Low-rank Matrix Approximations.
Jianwei Xiao, Ming Gu and Julien Langou.
24th IEEE International Conference on High Performance Computing, Data, and Analytics (HIPC), Jaipur, India, 2017.