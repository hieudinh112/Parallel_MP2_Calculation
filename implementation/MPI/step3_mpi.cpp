#include "reduction.h"
#include <cmath>
#include <iostream>
#include <mpi.h>
#include <chrono>

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);

    auto t_in = std::chrono::high_resolution_clock::now();
    
    // Inititalize the MPI World, get the size and the rank
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //// THESE PARAMETERS ARE PROBLEM SPECIFIC //////////
    // For testing purpose and reproducing the weak scaling analysis, the parameters ////////////
    // can be set as documented in the report //////////////////////////////////////////////////
    int NELECS = 30, NBASIS = 100, NAUX = 400;
    int N = 15; // For real application, we should have N = NELECS / size and the remainders 
    /////////////// should be distributed properly
    ////////////////////////////////////////////////////////////////////////////////
    // Small test: N: 15, NELECS: 30, NBASIS: 100, NAUX: 400.
    // Large test: N = 30, NELECS = 60, NBASIS = 200, NAUX = 800. 

    // Get the dimension of the problem
    if (rank == 0){
        std::cout << "========== Dimension of the problem ==========" << std::endl;\
        std::cout << "The size is: " << size << std::endl; 
        std::cout << "N: " << N << std::endl;
        std::cout << "NELECS: " << NELECS << std::endl;
        std::cout << "NBASIS: " << NBASIS << std::endl;
        std::cout << "NAUX: " << NAUX << std::endl;
    }

    // Initialize MP2 class and get the MP2 energy
    MP2 *correction = new MP2(NELECS, NBASIS, NAUX, N, size, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "========== Doing MP2 calculation ==========" << std::endl;
    }
    correction->get_energy();

    auto t_out = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = t_out - t_in;
    printf("Total calculation is done in %.2lfs\n", diff.count());
    std::cout << "====================================" << std::endl;

    delete correction;
    MPI_Finalize();

    return 0;
}
