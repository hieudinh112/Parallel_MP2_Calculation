#include <iostream>
#include <x86intrin.h>
#include <omp.h>

/** 
 * @brief This is the kernel to compute MP2 energy after obtaining rank-four tensor
 * @param NELECS number of electrons in the molecular system
 * @param NBASIS number of basis functions
 * @param ipjq pointer to the rank for matrix but index i and j have been fixed
 * @param MO_energy pointer ot the molecular orbital energy array
*/

// Cache locality can be an issue here.
double compute_MP2(const int NELECS,
                   const int NBASIS,
                   const int i,
                   const int j,
                   const double *ipjq,
                   const double *MO_energy)
{
    double E_MP2 = 0;
    int NVIR = NBASIS - NELECS;

    for (int p = 0; p < NVIR*NVIR; p++){
        // Requires further optimization
        int a = p/NVIR;
        int b = p%NVIR;
        double ibja = ipjq[b*NVIR + a];
        double iajb = ipjq[a*NVIR + b];
        double delta = MO_energy[i] + MO_energy[j] - MO_energy[a + NELECS] - MO_energy[b + NELECS];
        E_MP2 += 1/delta * (2 * iajb *iajb - iajb * ibja);
    }
    return E_MP2;
}