## Libraries requirement: 
gcc/12.1.0-fasrc01 OpenBLAS/0.3.7-fasrc01 openmpi/4.1.3-fasrc01

## Running the code
To execute the MPI for MP2 energy calcualtion kernel, run the following command
`sbatch submit_job.sh`

The step3_mpi can also be compiled by running
`make`
and the MPI can be run on a single node by mapping to different cores by running
`srun -n 4 -c 1 ./step3_mpi`

In here we also have a sample run 
- Small test: MPI_mp2_221893

By inspecting the *.out file, we can tell whether there's any work imbalance between MPI ranks.
The runtime can then be extracted, properly normalized, and used for weak scaling analysis
as being documented in the report. 