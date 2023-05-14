#include <iostream>
#include <chrono>
#include <fstream>
#include <omp.h>
#include <cblas.h>

int main (){

    // std::string MO_file = "./../../test/lih/lih_sto3g_MOenergy.out";
    // std::string coeff_file = "./../../test/lih/lih_sto3g_C.out";
    // std::string aux_file = "./../../test/lih/lih_sto3g_eri.out";
    // std::string overlap_file = "./../../test/lih/lih_sto3g_overlap.out";


    std::string MO_file = "./../../test/c3h8/c3h8_ccpvtz_MOenergy.out";
    std::string coeff_file = "./../../test/c3h8/c3h8_ccpvtz_C.out";
    std::string aux_file = "./../../test/c3h8/c3h8_eri.out";
    std::string overlap_file = "./../../test/c3h8/c3h8_overlap.out";
    /// Initialize parameters ///////////// We can change between vector and dynamic array
    size_t NELECS, NBASIS, NAUX;
    double SCF_ENERGY;

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
    ////////// STEP 0a:  READ IN THE MO energy ///////////////
    auto t0a_in = std::chrono::steady_clock::now();
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
    auto t0a_out = std::chrono::steady_clock::now();
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
    // TODO: openmp parallelization
    for (size_t P = 0; P < NAUX; P ++){
        for (size_t mu = 0; mu < NBASIS; mu++){
            for (size_t q = 0; q < NBASIS; q ++){
                for (size_t nu = 0; nu < NBASIS; nu++){
                    uqP[P * NBASIS * NBASIS + mu * NBASIS + q] += uvP[P * NBASIS * NBASIS + mu * NBASIS + nu] * C[q*NBASIS + nu];
                }
            }
        }
    }
    delete[] uvP;

    auto t1_mid = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff1_mid = t1_mid - t1_in;
    printf("Computing step 1a (MO Transformation) done in %.2lfs\n", diff1_mid.count());
    std::cout << "====================================" << std::endl;


    double *pqP = new double[NBASIS * NBASIS* NAUX]();

    // TODO: openmp parallelization
    for (size_t P = 0; P < NAUX; P ++){
        for (size_t p = 0; p < NBASIS; p ++){
            for (size_t q = 0; q < NBASIS; q++){
                for (size_t mu = 0; mu < NBASIS; mu++){
                    pqP[P * NBASIS * NBASIS + p * NBASIS + q] += uqP[P * NBASIS * NBASIS + mu * NBASIS + q] * C[p*NBASIS + mu];
                }
            }
        }
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
    // TODO: openmp parallelization
    for (size_t P = 0; P < NAUX; P ++){
        for (size_t Q = 0; Q < NAUX; Q ++){
            for (size_t p = 0; p < NBASIS; p++){
                for (size_t q = 0; q < NBASIS; q ++){
                    B_pqP[P * NBASIS * NBASIS + p * NBASIS + q] += pqP[Q * NBASIS * NBASIS + p * NBASIS + q] * overlap[P * NAUX + Q];
                }
            }
        }
    }

    auto t2_out = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff2 = t2_out - t2_in;
    printf("Computing step2 (Orthogonalization) done in %.2lfs\n", diff2.count());
    std::cout << "====================================" << std::endl;

    delete[] overlap;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////// STEP 3: MP2 energy calculation //////////////////////////////////////////////
    double E_MP2 = 0;
    auto t3_in = std::chrono::steady_clock::now();
    for (size_t i = 0; i < NELECS; i++){
        for (size_t j = 0; j < NELECS; j++){
            for (size_t a = NELECS; a < NBASIS; a++){
                for (size_t b = NELECS; b < NBASIS; b++){
                    double iajb = 0;
                    double ibja = 0;
                    for (size_t P = 0; P < NAUX; P++){
                        iajb += B_pqP[P * NBASIS * NBASIS + i * NBASIS + a] * B_pqP[P * NBASIS * NBASIS + j * NBASIS + b];
                        ibja += B_pqP[P * NBASIS * NBASIS + i * NBASIS + b] * B_pqP[P * NBASIS * NBASIS + j * NBASIS + a];
                    }
                    double delta = MO_energy[i] + MO_energy[j] - MO_energy[a] - MO_energy[b];
                    E_MP2 += 1/delta * (2 * iajb *iajb - iajb * ibja);
                }
            }
        }
    }
    auto t3_out = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff3 = t3_out - t3_in;
    printf("Computing step 3 (MP2 energy) done in %.2lfs\n", diff3.count());
    std::cout << "====================================" << std::endl;

    delete[] MO_energy;
    delete[] B_pqP;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "The final MP2 energy is: " << E_MP2 << std::endl;
    // The final energy should be -0.0174665 for LiH;
    // The final energy should be -1.1651e+14 for C3H8;
}