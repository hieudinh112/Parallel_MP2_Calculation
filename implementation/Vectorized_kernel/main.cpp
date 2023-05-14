#include <iostream>
#include <chrono>
#include <cblas.h>
#include <fstream>
#include <cstdlib>

extern double compute_MP2(const int NELECS,
                          const int NBASIS,
                          const int i,
                          const int j,
                          const double *ipjq,
                          const double *MO_energy);

extern "C" double compute_MP2_avx(const int NELECS,
                                  const int NBASIS,
                                  const int i,
                                  const int j,
                                  const double *ipjq,
                                  const double *MO_energy);

/////////////////////////////////////////////////////////////////////////////////////////
// In this script, I am skipping the basis transformation and orthogonalization step  //
////////////////////////////////////////////////////////////////////////////////////////

int main(void){
    
    // std::string MO_file = "./../../test/c3h8/c3h8_ccpvtz_MOenergy.out";
    // std::string coeff_file = "./../../test/c3h8/c3h8_ccpvtz_C.out";
    // std::string aux_file = "./../../test/c3h8/c3h8_eri.out";
    // std::string overlap_file = "./../../test/c3h8/c3h8_overlap.out";

    // int NELECS, NBASIS, NAUX;
    // double SCF_ENERGY;
    
    // Preprocess data
//     std::ifstream MO_mat;
//     MO_mat.open(MO_file);

//     MO_mat >> NELECS;
//     MO_mat >> NBASIS;
//     MO_mat >> NAUX;
//     MO_mat >> SCF_ENERGY;
    
//     std::cout << "NELECS: " << NELECS << std::endl;
//     std::cout << "NBASIS: " << NBASIS << std::endl;
//     std::cout << "NAUX: " << NAUX << std::endl;

//     int NVIR = NBASIS - NELECS;

//     double *MO_energy = new double[NBASIS];
//     for (int i = 0; i < NBASIS; i++){
//         MO_mat >> MO_energy[i];
//     }
//     MO_mat.close();

// ///////////////////////////////////////////////////////////
//     double *uvP = new double[NBASIS * NBASIS* NAUX];
//     std::ifstream eri;
//     eri.open(aux_file);
    
//     // #pragma omp parallel for num_threads(4) collapse(3) // This one her is wrong
//     // Need to rewrite this to make this faster
//     for (int P = 0; P < NAUX; P ++){
//         for (int mu = 0; mu < NBASIS; mu++){
//             for (int nu = 0; nu < NBASIS; nu++){
//                 eri >> uvP[P * NBASIS * NBASIS + mu * NBASIS + nu];
//             }
//         }
//     }
//     eri.close();
//     // std::precision(14);
//     std::cout << uvP[0] << std::endl;
// //////////////////////////////////////////////////////////

// ///////////////////////////////////////////////////////
//     // A bit of data post processing here so that the dimension works out
//     double *pqP_new = new double[NBASIS * NBASIS * NAUX] ();
//     for (int P = 0; P < NAUX; P++){
//         for (int p = 0; p < NBASIS; p++){
//             for (int q = 0; q < NBASIS; q++){
//                 pqP_new[p*NAUX*NBASIS + q *NAUX + P] = uvP[P*NBASIS*NBASIS + p * NBASIS + q];
//             }
//         }
//     }
//     delete[] uvP;

//     double *B_iaP = new double[NVIR * NELECS * NAUX] ();
//     for (int P = 0; P < NAUX; P++){
//         for (int i = 0; i < NELECS; i ++){
//             for (int a = NELECS; a < NBASIS; a++){
//                 B_iaP[i * NVIR * NAUX + P * NVIR + a - NELECS] = pqP_new[i*NAUX*NBASIS + a *NAUX + P];
//             }
//         }
//     }
//     delete [] pqP_new;

    ///// RANDOM INITIALIZATION /////////////
    int NELECS = 50, NBASIS = 300, NAUX = 800; // THESE VARIABLES ARE SPECIFIC FOR MOLECULE 
    int NVIR = NBASIS - NELECS;
    /////////////////////////////////////////

    double *B_iaP = new double[NVIR * NELECS * NAUX]();
    double *MO_energy = new double[NBASIS];

    for (int i = 0; i < NBASIS; i++){
        MO_energy[i] = i * 0.5; // just let the MO_energy be the same as i
    }
    
    // double *B_iaP = new double[NVIR * NELECS * NAUX] ();
    for (int P = 0; P < NAUX; P++){
        for (int i = 0; i < NELECS; i ++){
            for (int a = NELECS; a < NBASIS; a++){
                B_iaP[i * NVIR * NAUX + P * NVIR + a - NELECS] = (double) rand()/RAND_MAX;
            }
        }
    }
    ///////////////////////////////////////////////////////////////////////

    ////////////////// BENCHMARK INSIDE THIS REGION //////////////////////
    // auto t3_in = std::chrono::steady_clock::now();
    auto t3_in = std::chrono::steady_clock::now();
    double E_MP2 = 0;
    for (int i = 0; i < NELECS; i++){
        for (int j = 0; j < NELECS; j++){
            double *ipjq = new double[NVIR * NVIR] ();
            // contract over P
            cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, NVIR, NVIR, NAUX, 1.0, B_iaP + i*NVIR*NAUX, NVIR, B_iaP+j*NVIR*NAUX, NVIR, 1.0, ipjq, NVIR);
            // Compare the vectorize kernel here
            E_MP2 += compute_MP2(NELECS, NBASIS, i, j, ipjq, MO_energy);
            delete [] ipjq;
        }
    }
    auto t3_mid = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_mid = t3_mid - t3_in;
    printf("Computing step 3 (MP2 energy) done in %.2lfs\n", time_mid.count());


    double E_MP2_avx2 = 0;
    for (int i = 0; i < NELECS; i++){
        for (int j = 0; j < NELECS; j++){
            double *ipjq = new double[NVIR * NVIR] ();
            // contract over P
            cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, NVIR, NVIR, NAUX, 1.0, B_iaP + i*NVIR*NAUX, NVIR, B_iaP+j*NVIR*NAUX, NVIR, 1.0, ipjq, NVIR);
            // Compare the vectorize kernel here
            E_MP2_avx2 += compute_MP2_avx(NELECS, NBASIS, i, j, ipjq, MO_energy);
            delete [] ipjq;
        }
    }

    auto t3_out = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_diff = t3_out - t3_mid;
    printf("Computing step 3 (MP2 energy) (vectorized) done in %.2lfs\n", time_diff.count());
    std::cout << "====================================" << std::endl;

    delete [] B_iaP;
    delete [] MO_energy;


    // This part is for checking the consistency
    double diff = std::abs(E_MP2 - E_MP2_avx2);
    // std::cout << "E_MP2 = " << E_MP2 << ", E_MP2_AVX = " << E_MP2_avx2 << std::endl;
    std::cout << "====== Exact energy difference is: " << diff << " ======" << std::endl;

    return 0;
}