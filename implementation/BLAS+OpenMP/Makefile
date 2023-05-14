# Make File to compile and run the mp2 calculation
CXX = g++
CXXFLAGS = -Wall -std=c++11 -fopenmp

.PHONY: run clean

ri_mp2_blas_ver2: ri_mp2_blas_ver2.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -lopenblas

ri_mp2: ri_mp2.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

ri_mp2_blas: ri_mp2_blas.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -lopenblas

ri_mp2_blas_ver2_omp: ri_mp2_blas_ver2_omp.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -lopenblas
clean:
	rm -f ri_mp2 ri_mp2_blas ri_mp2_blas_ver2 ri_mp2_blas_ver2_omp

run: ri_mp2
	 ./ri_mp2
