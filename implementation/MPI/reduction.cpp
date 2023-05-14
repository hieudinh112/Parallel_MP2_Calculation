#include "reduction.h"
#include <cmath>
#include <cblas.h>
#include <chrono>
#include <cstdlib> 


MP2::MP2(const int NELECS,
         const int NBASIS,
         const int NAUX,
         const int N, // Number of occupied orbitals per rank
         const int proc_num, // number of processors in MPI WORLD
         const MPI_Comm comm)
    : NELECS_(NELECS), NBASIS_(NBASIS), NAUX_(NAUX), NVIR_(NBASIS - NELECS), proc_num_(proc_num), N_(N) //, MO_energy_(MO_energy)
{
    // Initialize the Cartesian topology
    period_[0] = true;
    int dims[1] = {proc_num_};
    MPI_Cart_create(comm, 1, dims, period_, true, &cart_comm_); // int MPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[], const int periods[], int reorder, MPI_Comm * comm_cart)

    // Determine the neighbor in the Cartesian topology
    MPI_Cart_shift(cart_comm_, 0, 1, &neighbor_[Left], &neighbor_[Right]); // int MPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest)

    // Determine the rank of the current proccesor in the Cartesian topology
    MPI_Comm_rank(cart_comm_, &cart_rank_);

    idx_A = cart_rank_;
    idx_B = cart_rank_;

    // Initialize the matrix here
    A_ = new double[N_*NVIR_*NAUX_];
    B_ = new double[N_*NVIR_*NAUX_];
    Brecv_ = new double[N_*NVIR_*NAUX_] (); // Initialize to 0
    MO_energy_ = new double[NBASIS_];

    initialize_(); // Put the data into A and B 
}

void MP2::initialize_(){
    //random initialization for matrices A_ and B_
    for (int i = 0; i < N_; i++){
        for (int P = 0; P < NAUX_; P++){
            for (int a = 0; a < NVIR_; a++){
                double element = 1 * 1e-5 + cart_rank_ * 1e-4 * 5 + (i +P + a) * 1e-6;
                A_[i * NVIR_ * NAUX_ + P * NVIR_ + a] = element;
                B_[i * NVIR_ * NAUX_ + P * NVIR_ + a] = element;
            }
        }
    }

    for (int i = 0; i < NBASIS_; i++){
        MO_energy_[i] = i * 0.5; // just let the MO_energy be the same as i
    }
    MP2_energy = 0;
    MP2_partial = 0;
}

// This function is used to make sure that the correct MO energy is used (only need to update the idx_B for memory at B)
void MP2::update_index(){
    idx_B -= 1;
    if (idx_B < 0){
        idx_B += proc_num_;
    }
}

// This is the compute kernel to calculate the MP2 energy (last step in our pipeline)
void MP2::reduction_kernel(){
    // Each rank only has acess to n elements of the 
    for (int i = 0; i < N_; i++){
        for (int j = 0; j < N_; j++){
            double *ipjq = new double[NVIR_ * NVIR_] ();
            
            cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, NVIR_, NVIR_, NAUX_, 1.0, A_ + i*NVIR_*NAUX_, NVIR_, B_ + j*NVIR_*NAUX_, NVIR_, 1.0, ipjq, NVIR_); // A and B holds some portions of rank three tensor
            
            for (int a = NELECS_; a < NBASIS_; a++){
                for (int b = NELECS_; b < NBASIS_; b++){
                    double ibja = ipjq[(b - NELECS_)*NVIR_ + (a - NELECS_)];
                    double iajb = ipjq[(a - NELECS_)*NVIR_ + (b - NELECS_)];
                    double delta = MO_energy_[i + idx_A*N_] + MO_energy_[j + idx_B*N_] - MO_energy_[a] - MO_energy_[b]; // Need to change the rank ID as well (must be able to figure this out)
                    MP2_partial += 1/delta * (2 * iajb *iajb - iajb * ibja);
                }
            }
            delete [] ipjq;
        }
    }
}

// This function contains the implementation of compute/transfer overlap procedure
// Inspired by the amazing Cannon algorithm for matrix matrix multiplication
void MP2::get_energy(){
    MPI_Request request[2];
    for (int k = 1; k < proc_num_; ++k){
        MPI_Isend(B_, NAUX_*NVIR_*N_, MPI_DOUBLE, neighbor_[Right], 123, cart_comm_, &request[0]); 
        MPI_Irecv(Brecv_, NAUX_*NVIR_*N_, MPI_DOUBLE, neighbor_[Left], 123, cart_comm_, &request[1]);
        reduction_kernel();
        update_index();
        MPI_Waitall(2, request, MPI_STATUS_IGNORE);
        std::swap(B_, Brecv_);
    }
    // The final calculation is done outside of the loop
    reduction_kernel();
    MPI_Reduce(&MP2_partial, &MP2_energy, 1, MPI_DOUBLE, MPI_SUM, 0, cart_comm_);
}

// Destructor
MP2::~MP2(){
    // Free the allocated energy after MP2 calculation
    delete[] A_;
    delete[] B_;
    delete[] Brecv_;
    delete[] MO_energy_;
    // Free the MPI Cartesian communicator
    MPI_Comm_free(&cart_comm_);
    }