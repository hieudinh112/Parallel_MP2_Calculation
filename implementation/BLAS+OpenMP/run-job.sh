#!/usr/bin/env bash
#SBATCH --job-name=omp_test_%j
#SBATCH --output=mp2_%j.out
#SBATCH --error=jmp2_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:30:00
# set the number of OpenMP threads per MPI process
# (more on this in class later)


make clean
make ri_mp2_blas_ver2_omp

for N in 1 2 4 8 16 32
do 
    echo -e "\n*************************** \n Run with ${N} threads: \n"

    OMP_PROC_BIND='true' OMP_NUM_THREADS=${N} srun ./ri_mp2_blas_ver2_omp "$@"           
done
# return the exit code of srun above
exit $?
