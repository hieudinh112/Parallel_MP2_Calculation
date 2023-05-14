#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <omp.h>
#include <iomanip>
#include <cmath>
#include <cstdio>

// all calulations done here is restricted calculation

// void read_MO_energy(std::string pathname, energy_array);

// void read_2e3c_integral(std::string pathname, std::vector<double>);

void read_MO_energy(std::string pathname, double &energy[], int &NELECS, int &NBASIS){
    
    std::ifstream input(pathname);
    
    int NELECS, NBASIS, NAUX; 
    double SCF_ENERGY;
    
    input >> NELECS;
    input >> NBASIS;
    input >> NAUX;
    input >> SCF_ENERGY;

    for (int i = 0; i < NBASIS; ++i){
        input >> energy[i];
    }
    input.close();
}

// void read_2e3c_ints(std::string pathname, double &B_mat[]){

// }

int get_MP2(double B_mat[], double MO_energy[]){
    
    int NBASIS = 7;
    int NELECS = 5;
    int NAUX = 4;

    double mp2_energy = 0;

    auto t0 = std::chrono::steady_clock::now();

    // run time is O(Nvir^2 Nocc^2 Naux)
    for (int i = 0; i < NELECS; ++i){
        for (int j = 0; j < NBASIS; ++j){
            for (int a = NELECS; a < NBASIS; ++a){
                for (int b = NELECS; b < NBASIS; ++b){
                    for (int P = 0; P < NAUX; ++P){
                        iajb = B_mat[i][a][P] * B_mat[j][b][P]; // we can potentially exchange loop to optimize this multiplication
                        ibja = B_mat[i][b][P] * B_mat[j][a][P];
                    }
                    delta = MO_energy[i] + MO_energy[j] - MO_energy[a] - MO_energy[b];
                    mp2 += iajb * (2 * iajb - ibja) / delta;
                }
            }
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = t1 - t0;
    printf("Computing done in %.2lfs\n for get_MP2 function", diff.count());
    std::cout << "The final mp2 energy calculation is: " << mp2 << std::endl;
    return 0;
}

int main(){
    
    // std::string energy_path = "";
    // std::string 2e3c_path = "";

    int NELECS, NBASIS, NAUX; 
    double SCF_ENERGY;

    double MO_energy[];
    read_MO_energy(energy_path, MO_energy, NELECS, NBASIS, NAUX, SCF_ENERGY);

    std::cout << "NELECS" << NELECS << std::endl;
    std::cout << "NBASIS" << NBASIS << std::endl;
    std::cout << "NAUX" << NAUX << std::endl;
    std::cout << "SCF_ENERGY" << SCF_ENERGY << std::endl;

    for (int i =0; i < sizeof(MO_energy) / sizeof(MO_energy[0]); ++i){
        std::cout << "MO_energy(i)" << MO_energy[i] << std::endl;
    }

    return 0;

}
