#include <iostream>
#include <chrono>
#include <fstream>
#include <omp.h>
#include <mpi.h>
#include <cblas.h>

int main(int argc, char *argv[]){

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    std::string MO_file = "/n/academic_homes/g113711/u449542g113711/team10/test/lih/lih_sto3g_MOenergy.out";
    std::string coeff_file = "/n/academic_homes/g113711/u449542g113711/team10/test/lih/lih_sto3g_C.out";
    std::string aux_file = "/n/academic_homes/g113711/u449542g113711/team10/test/lih/lih_sto3g_eri.out";
    std::string overlap_file = "/n/academic_homes/g113711/u449542g113711/team10/test/lih/lih_sto3g_overlap.out";
    
    size_t NELECS, NBASIS, NAUX;
    double SCF_ENERGY;


    // Step 0: Read in the MO energy
    MO_mat.open(MO_file);
    MO_mat >> NELECS;
    MO_mat >> NBASIS;
    MO_mat >> NAUX;
    MO_mat >> SCF_ENERGY;
    
    std::cout << "NELECS: " << NELECS << std::endl;
    std::cout << "NBASIS: " << NBASIS << std::endl;
    std::cout << "NAUX: " << NAUX << std::endl;

    double *MO_energy = new double[NBASIS];
    for (size_t i = 0; i < NBASIS; i++){
        MO_mat >> MO_energy[i];
    }
    MO_mat.close();

    // Step 0: Reading in the coeff matrix (all ranks read)
    double *C = new double[NBASIS * NBASIS];
    std::ifstream coeff_mat;
    coeff_mat.open(coeff_file);
    
    for (size_t p = 0; p < NBASIS; p++){
        for (size_t nu = 0; nu < NBASIS; nu++){
            coeff_mat >> C[p*NBASIS + nu]; //(We are reading in C_pnu)
        }
    }
    coeff_mat.close();

    // Step 0: Read in the overlap matrix (all ranks read)
    std::ifstream overlap_mat;
    overlap_mat.open(overlap_file);

    double *overlap = new double[NAUX * NAUX];
    for (size_t P = 0; P < NAUX; P++){
        for (size_t Q = 0; Q < NAUX; Q++){
            overlap_mat >> overlap[P* NAUX + Q];
        }
    }
    overlap_mat.close();

    // Step 0: Read in rank 3 tensor from the root rank and then scatter it to other ranks

    if (rank == 0){
        double *uvP = new double[NBASIS * NBASIS* NAUX];
        std::ifstream eri;
        eri.open(aux_file);
        
        // #pragma omp parallel for num_threads(4) collapse(3) // This one her is wrong
        // Need to rewrite this to make this faster
        for (size_t P = 0; P < NAUX; P ++){
            for (size_t mu = 0; mu < NBASIS; mu++){
                for (size_t nu = 0; nu < NBASIS; nu++){
                    eri >> uvP[P * NBASIS * NBASIS + mu * NBASIS + nu];
                }
            }
        }
        eri.close();
    } else {

    }



    // Step 1:


    // Step 2:



    // Step 3:



    MPI_Finalize();
}