#ifndef REDUCTION_H
#define REDUCTION_H

#include <iostream>
#include <mpi.h>

class MP2
{

public:
    /**
     * @brief Constructor of MP2 class
     * @param NELECS Number of electrons
     * @param NBASIS Number of basis functions
     * @param NAUX Number of auxiliary functions
     * @param N Number of occupied orbitals per rank (This number needs to be fixed to do weak scaling analysis)
     * @param MO_energy pointer to the array of MO orbitals
     * @param proc_num Number of MPI processes in the calculation
     * @param comm Global MPI communicator
     */
    MP2(const int NELECS,
        const int NBASIS,
        const int NAUX,
        const int N,
        // const double *MO_energy,
        const int proc_num,
        const MPI_Comm comm = MPI_COMM_WORLD);

    void reduction_kernel();
    void update_index();
    void get_energy();
    void initialize_();
    ~MP2();

private:
    // Dimensions relevant to the MP2 class (specific molecule)
    int NELECS_, NBASIS_, NAUX_, NVIR_;

    int idx_A, idx_B; // These indices are used to index into the MO_energy
    // MPI rank neighbor enumeration
    enum Neighbor {Left, Right, NNeighbors};
    // MPI rank in Cartesian topology
    int cart_rank_;
    // Others
    int neighbor_[NNeighbors], period_[1], coords_;
    // Cartesian MPI communicator
    MPI_Comm cart_comm_;

    const int proc_num_; // Number of processors
    const int N_; // Number of occupied orbitals in this rank

    // Store the rank 3 matrix
    double *A_;
    double *B_;
    double *Brecv_;
    double *MO_energy_;

    // MP2 energy
    double MP2_partial; // Partial MP2 energy in each rank
    double MP2_energy;
};

#endif // REDUCTION_H


