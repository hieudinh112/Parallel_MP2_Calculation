#include <omp.h>
#include <iostream>
#include <fstream>
#include <strings>
using namespace std;



// read in both text files
//https://stackoverflow.com/questions/36708370/reading-from-txt-file-into-two-dimensional-array-in-c
void read_in_B(void){
	ifstream in;
	in.open("B_matrix.txt")
	float element;


	if (in.is_open()){
//ToDo: figure out correct chunk size
#pragma omp parallel for collapse(3) schedule(static, 1)
		for(int i=0, i<N, i++){
			for(int j =0, j<M, j++){
				for(int k=0, k<P, k++){
					in >> B_mat[i][j][k];
				}
			}
		}
	}

}

void read_in_energies(void){
	ifstream in;
	in.open("energies.txt")
	float element;


	if (in.is_open()){
		for(int i=0, i<N_basis, i++){
				in >> energies[i];
		}
	}

}






void eri_tensor(void){
//ToDo: figure out correct chunk size, can we use reduction?
// reduction resource: https://stackoverflow.com/questions/20413995/reducing-on-array-in-openmp
#pragma omp parallel for collapse(5) schedule(static, 1) 
	for(int i=0, i<N, i++){
		for(int j=0, j<N, j++){
			for(int a=N, a < N_basis, a++){
				for(int b=N, b< N_basis, b++){
					for(int p=0, p<P, p++){
						eri[i][a][j][b] += B[i][a][p]*B[j][b][p];
					}
				}
			}
		}
	}
}

//compute energy
float main(void){

	// change number of basis functions for desired application (occ,vir,aux)
	int N = 8;
	int M = 8;
	int P = 80
	int N_basis = N+M;

	// initialize necessary parameters
	//ToDO: figure out the way QChem outputs the B tensor and fix the read in skeleton appropriately
	float B_mat[N][M][P];
	//ToDO: figure out the way QChem outputs the e array and fix the read in skeleton appropriately
	float energies[N_basis];

	read_in_B();
	read_in_energies();

	//build four index tensor
	float eri[N][M][N][M] = {0.0};
	eri_tensor();

	//initialize energy
	float E_MP2 = 0;

//ToDo: figure out correct chunk size
#pragma omp parallel for collapse(4) schedule(static, 1) reduction(+:E_MP2)
	for(int i =0, i<N, i++){
		for(int j =0, j<N, j++){
			for(int a=N, a<N_basis, a++){
				for(int b=N, b<N_basis,  b++){
					float delta = energies[i]+energies[j]-energies[a]-energies[b];
					E_MP2 += eri[i][a][j][b]*(2*eri[i][a][j][b]-eri[i][b][j][a])/delta;
				}
			}
		}
	}
	return E_MP2
}



