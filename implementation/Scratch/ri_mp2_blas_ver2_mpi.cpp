#include <iostream>
#include <chrono>
#include <fstream>
#include <omp.h>
#include <mpi.h>
#include <cblas.h>


int main (int argc, char **argv){

    std::string MO_file = "/n/academic_homes/g113711/u449542g113711/team10/test/lih/lih_sto3g_MOenergy.out";
    std::string coeff_file = "/n/academic_homes/g113711/u449542g113711/team10/test/lih/lih_sto3g_C.out";
    std::string aux_file = "/n/academic_homes/g113711/u449542g113711/team10/test/lih/lih_sto3g_eri.out";
    std::string overlap_file = "/n/academic_homes/g113711/u449542g113711/team10/test/lih/lih_sto3g_overlap.out";

    // std::string MO_file = "/n/academic_homes/g113711/u449542g113711/team10/test/c3h8/c3h8_ccpvtz_MOenergy.out";
    // std::string coeff_file = "/n/academic_homes/g113711/u449542g113711/team10/test/c3h8/c3h8_ccpvtz_C.out";
    // std::string aux_file = "/n/academic_homes/g113711/u449542g113711/team10/test/c3h8/c3h8_eri.out";
    // std::string overlap_file = "/n/academic_homes/g113711/u449542g113711/team10/test/c3h8/c3h8_overlap.out";

    /// Initialize parameters ///////////// We can change between vector and dynamic array
    size_t NELECS, NBASIS, NAUX;
    double SCF_ENERGY;

    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    ////////// STEP 0a:  READ IN THE MO energy ///////////////
    auto t0a_in = std::chrono::high_resolution_clock::now();
    std::ifstream MO_mat;
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
    auto t0a_out = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff0a = t0a_out - t0a_in;
    printf("Reading in MO energy done in %.2lfs\n", diff0a.count());
    std::cout << "====================================" << std::endl;
    ////////////////////////////////////////////////////////////////
    
    //////////////// STEP 0b: READ IN OVERLAP MAT ///////////////////////////
    auto t0b_in = std::chrono::steady_clock::now();
    std::ifstream overlap_mat;
    overlap_mat.open(overlap_file);

    double *overlap = new double[NAUX * NAUX];
    for (size_t P = 0; P < NAUX; P++){
        for (size_t Q = 0; Q < NAUX; Q++){
            overlap_mat >> overlap[P* NAUX + Q];
        }
    }
    overlap_mat.close();

    auto t0b_out = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff0b = t0b_out - t0b_in;
    printf("Reading in overlap matrix done in %.2lfs\n", diff0b.count());
    std::cout << "====================================" << std::endl;
    // std::cout << overlap[0] << std::endl;
    /////////////////////////////////////////////////////////////////

    /////////// STEP 0c: READ IN COEFFICIENT MATRIX ////////////
    auto t0c_in = std::chrono::steady_clock::now();
    double *C = new double[NBASIS * NBASIS];
    std::ifstream coeff_mat;
    coeff_mat.open(coeff_file);
    
    for (size_t p = 0; p < NBASIS; p++){
        for (size_t nu = 0; nu < NBASIS; nu++){
            coeff_mat >> C[p*NBASIS + nu]; //(We are reading in C_pnu)
        }
    }

    coeff_mat.close();
    auto t0c_out = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff0c = t0c_out - t0c_in;
    printf("Reading in coefficient matrix done in %.2lfs\n", diff0c.count());
    std::cout << "====================================" << std::endl;
    ///////////////////////////////////////////////////////////

    ///////////// STEP 0d: READ IN THE (uv|P) MATRIX //////////
    auto t0d_in = std::chrono::steady_clock::now();
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

    auto t0d_out = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff0d = t0d_out - t0d_in;
    printf("Reading in rank-3 tensor done in %.2lfs\n", diff0d.count());
    std::cout << "====================================" << std::endl;
    /////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////// STEP 1: MO TRANSFORMATION //////////////////////////////////////////
    double *uqP = new double[NBASIS * NBASIS* NAUX] ();

    auto t1_in = std::chrono::steady_clock::now();
    // TODO: openmp/MPI parallelizatio
     MPI_Init(&argc,&argv);

     //ranks and size
     int rank, size;
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &size);
     printf("size %d. \n", size);

    //figure out how many elements in to send to each processor (like lab4)
    int n_elements = NAUX/size;
    int remainder = NAUX%size;

    int first = n_elements*rank;
    int last = n_elements*(rank+1)-1;

    //distribute remainder
    if (rank%size<remainder){
          first += rank;
          last += rank+1;
      } else{
          first += remainder;
          last += remainder; 
      }
    
    int num_elements = last-first+1;

    //make array so we know how many elements in each processor
    int recvcounts[size];
    MPI_Allgather(&num_elements,1,MPI_INT, &recvcounts,1, MPI_INT, MPI_COMM_WORLD);
    //make an array to receive displacements
    int *displacements = new int[size];
    int cbytes = 8*NBASIS*NBASIS*(num_elements);
    MPI_Exscan(&cbytes, &displacements, NBASIS*NBASIS*(num_elements), MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    printf("num elements %d. \n", num_elements);
    //create intermediate matrices
    //double *uvP_intermediate = new double[NBASIS * NBASIS* (last-first)] ();
    double *uqP_intermediate = new double[NBASIS * NBASIS* (last-first)] ();

    // int left, right;
    // right = (rank + 1) % size;
    // left = rank - 1;
    // if (left < 0)
    //     left = size - 1;
    // //use nonblocking communication to send each chunk
    // MPI_Request request[size];

    // MPI_Isend(&uvP+first*NBASIS*NBASIS, (last-first)*NBASIS*NBASIS, MPI_DOUBLE,right, 123,MPI_COMM_WORLD,&request[0]);
    // MPI_Irecv(&uvP_intermediate, (last-first)*NBASIS*NBASIS, MPI_DOUBLE,left, 123,MPI_COMM_WORLD,&request[1]);
    // MPI_Waitall(2, request, MPI_STATUSES_IGNORE);

    //perform computation 
    for (int P = first; P < last; P ++){
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, NBASIS, NBASIS, NBASIS, 1.0, uvP + P *NBASIS *NBASIS, NBASIS, C, NBASIS, 1.0, uqP_intermediate+ (P-first) *NBASIS *NBASIS, NBASIS);
    }
    // for (size_t P = 0; P < NAUX; P ++){
    //     cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, NBASIS, NBASIS, NBASIS, 1.0, uvP + P *NBASIS *NBASIS, NBASIS, C, NBASIS, 1.0, uqP+ P *NBASIS *NBASIS, NBASIS);
    // }
    //send the results back to master thread
    
    //gather results
    MPI_Gatherv(&uqP_intermediate, (last-first)*NBASIS*NBASIS, MPI_DOUBLE, &uqP, recvcounts,displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Finalize();
    delete[] uvP;
    auto t1_mid = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff1_mid = t1_mid - t1_in;
    printf("Computing step 1a (MO Transformation) done in %.2lfs\n", diff1_mid.count());
    std::cout << "====================================" << std::endl;


    double *pqP = new double[NBASIS * NBASIS* NAUX]();

    // TODO: openmp/MPI parallelization
    for (size_t P = 0; P < NAUX; P ++){
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NBASIS, NBASIS, NBASIS, 1.0, C, NBASIS, uqP + P *NBASIS *NBASIS, NBASIS, 1.0, pqP+ P *NBASIS *NBASIS, NBASIS);
    }
    auto t1_out = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff1 = t1_out - t1_mid;
    printf("Computing step 1b (MO Transformation) done in %.2lfs\n", diff1.count());
    std::cout << "====================================" << std::endl;

    delete[] uqP;
    delete[] C;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////// STEP2: ORTHOGONALIZATION /////////////////////////////////////////////////////
    double *B_pqP = new double[NBASIS * NBASIS * NAUX] ();

    auto t2_in = std::chrono::steady_clock::now();
    // Reshaping for pqP //  
    double *pqP_new = new double[NBASIS * NBASIS * NAUX] ();

    for (size_t P = 0; P < NAUX; P++){
        for (size_t p = 0; p < NBASIS; p++){
            for (size_t q = 0; q < NBASIS; q++){
                pqP_new[p*NAUX*NBASIS + q *NAUX + P] = pqP[P*NBASIS*NBASIS + p * NBASIS + q];
                // pqP_new[p*NAUX*NBASIS + P*NBASIS + q] = pqP[P*NBASIS*NBASIS + p * NBASIS + q];
            }
        }
    }
    delete[] pqP;

    // TODO: openmp/MPI parallelization (maybe)
    // overlap is a symmetric matrix so row major order or column major order do not matter
    for (size_t p = 0; p < NBASIS; p++){
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NBASIS, NAUX, NAUX, 1.0, pqP_new+p*NAUX*NBASIS, NAUX, overlap, NAUX, 1.0, B_pqP+p*NAUX*NBASIS, NAUX);
        // cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, NBASIS, NAUX, NAUX, 1.0, pqP_new+p*NAUX*NBASIS, NBASIS, overlap, NAUX, 1.0, B_pqP+p*NAUX*NBASIS, NAUX);
    }
    auto t2_out = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff2 = t2_out - t2_in;
    printf("Computing step2 (Orthogonalization) done in %.2lfs\n", diff2.count());
    std::cout << "====================================" << std::endl;

    delete[] pqP_new;
    delete[] overlap;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////// STEP 3: MP2 energy calculation //////////////////////////////////////////////
    double E_MP2 = 0;
    auto t3_in = std::chrono::steady_clock::now();

    ////// Reshaping of the B matrix to keep only ov block (i is the slowest changing index here) //////
    double * B_iaP = new double[(NBASIS - NELECS) * NELECS *NAUX] ();
    //openmp parallelization here?
    for (size_t P = 0; P < NAUX; P++){
        for (size_t i = 0; i < NELECS; i ++){
            for (size_t a = NELECS; a < NBASIS; a++){
                B_iaP[i * (NBASIS - NELECS) * NAUX + P * (NBASIS - NELECS) + a - NELECS] = B_pqP[i*NAUX*NBASIS + a *NAUX + P];
            }
        }
    }
    delete [] B_pqP;
    // auto t3_in = std::chrono::steady_clock::now();
    for (size_t i = 0; i < NELECS; i++){
        for (size_t j = 0; j < NELECS; j++){
            int NVIR = NBASIS - NELECS;
            double *ipjq = new double[NVIR * NVIR] ();
            // contract over P
            cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, NVIR, NVIR, NAUX, 1.0, B_iaP + i*NVIR*NAUX, NVIR, B_iaP+j*NVIR*NAUX, NVIR, 1.0, ipjq, NVIR);
            // contract over a and b
            // Vectorization can be done here since this is matrix matrix element wise multiplication
            for (size_t a = NELECS; a < NBASIS; a++){
                for (size_t b = NELECS; b < NBASIS; b++){
                    double ibja = ipjq[(b - NELECS)*NVIR + (a - NELECS)];
                    double iajb = ipjq[(a - NELECS)*NVIR + (b - NELECS)];
                    double delta = MO_energy[i] + MO_energy[j] - MO_energy[a] - MO_energy[b];
                    E_MP2 += 1/delta * (2 * iajb *iajb - iajb * ibja);
                }
            }
            delete [] ipjq;
        }
    }
    auto t3_out = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff3 = t3_out - t3_in;
    printf("Computing step 3 (MP2 energy) done in %.2lfs\n", diff3.count());
    std::cout << "====================================" << std::endl;

    delete[] MO_energy;
    delete[] B_iaP;
    // delete[] B_pqP;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "The final MP2 energy is: " << E_MP2 << std::endl;
    // The final energy should be -0.0174665 for LiH;
    // The final energy should be -1.1651e+14 for C3H8;
}