#include <iostream>
#include <chrono>
#include <fstream>
#include <omp.h>
#include <cblas.h>

int test_case(){

    double *A1 = new double[4 * 4]{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    double *C1 = new double[4*4]();
    std::cout << *(A1 + 2) << std::endl;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 2,2,2, 1.0, A1 + 2,2, A1+2, 2,1.0, C1, 2);
    // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 2,2,2, 1.0, A1 + 2,2, A1+2, 2,1.0, C1, 2);
    // cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasTrans, 4,4,4, 1.0, A1, 4, A1, 4, 1.0, C1, 4); // This is A1 * A1.T
    // cblas_dgemm(CblasColMajor, CblasNoTrans,CblasTrans, 4,4,4, 1.0, A1, 4, A1, 4, 1.0, C1, 4);
    // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,);
    for (size_t i = 0; i < 4*4; i++){
        std::cout << "C1[" << i << "]: " << C1[i] << std::endl;
    }
    delete[] A1;
    delete[] C1;
    return 0;
}

int main(){
    return 0;
}