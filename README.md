# psycloned_nemo

This repo simply contains the PSyclone-processed NEMO source files using the Met Office branch NEMO_4.0_mirror_SI3_GPU_13325

This branch contains PSyclone-processed code with OpenMP applied to loops over vertical levels (loop index 'jk' in NEMO). 
No dependence analysis was used in the processing step, so manual tweaks are needed to isolate race conditions.

The preprocessed source was generated using ./makenemo with MPI and SI3 keys added. The source was then processed 
using PSyclone branch nemo_mpi_tryout, and an OpenMP transformation script taken from the PSyclone tutorial:
https://github.com/stfc/PSyclone/blob/master/tutorial/practicals/nemo/3_nemo_openmp/solutions/general_parallel_region_omp_trans.py

Note that YOUR_ARCH_FILE must contain appropriate flags to compile with OpenMP enabled (-qopenmp for Intel). 
