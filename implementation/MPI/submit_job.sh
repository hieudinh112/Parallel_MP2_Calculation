#!/usr/bin/env bash
#SBATCH --job-name=mp2_energy
#SBATCH --output=MPI_mp2_%j.out
#SBATCH --error=MPI_mp2_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00

module purge
module load gcc/12.1.0-fasrc01 openmpi/4.1.3-fasrc01 OpenBLAS/0.3.7-fasrc01
make
## The number of rank here can be changed
for p in 1 2 4 8 16 32
do
    echo "========== Running MPI with $p ranks =========="
    srun -n $p  -c 1 ./step3_mpi
done 
make clean
