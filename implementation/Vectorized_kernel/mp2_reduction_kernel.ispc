// This file contains the ISPC implementation for the compute kernel

/** 
 * @brief This is the kernel to compute MP2 energy after obtaining rank-four tensor
 * @param NELECS number of electrons in the molecular system
 * @param NBASIS number of basis functions
 * @param ipjq pointer to the rank for matrix but index i and j have been fixed
 * @param MO_energy pointer ot the molecular orbital energy array
*/
export uniform double compute_MP2_avx(const uniform int NELECS,
                              const uniform int NBASIS,
                              const uniform int i,
                              const uniform int j,
                              const double *uniform ipjq,
                              const double *uniform MO_energy)
{
    varying double E_MP2 = 0;
    uniform int NVIR = NBASIS - NELECS;

    foreach (p = 0 ... NVIR*NVIR){
        varying int a = p/NVIR;
        varying int b = p%NVIR;
        varying double ibja = ipjq[b*NVIR + a];
        varying double iajb = ipjq[a*NVIR + b];
        varying double delta = MO_energy[i] + MO_energy[j] - MO_energy[a + NELECS] - MO_energy[b + NELECS];
        E_MP2 += 1/delta * (2 * iajb *iajb - iajb * ibja);
    }
    return reduce_add(E_MP2);
}