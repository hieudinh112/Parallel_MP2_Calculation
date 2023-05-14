## Libraries Requirement
gcc/12.1.0-fasrc01 OpenBLAS/0.3.7-fasrc01

## Files
- ri_mp2.cpp: Naive Implementation
- ri_mp2_blas.cpp: BLAS Implemetation (Version 1)
- ri_mp2_blas_ver2.cpp: BLAS Implementation (Version 2)
- ri_mp2_blas_ver2_omp.cpp: OpenMP Parallelization

## Run
To run the file, type in the following command
- `make file_name_no_extension`
- `./ri_mp2*`

The OpenMP parallelization script can be run as follow 
`sbatch run-job.sh`
and the runtime can be extracted to make the strong scaling plot
