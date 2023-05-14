# This is to compile the code for MPI calculation

CXX = mpic++
CXXFLAGS = -Wall -std=c++11 -fopenmp -Wextra -Wpedantic -lopenblas

.PHONY: run clean

SRC = step3_mpi.cpp reduction.cpp
DEP = reduction.h

all: step3_mpi

step3_mpi: $(SRC) $(DEP)
	$(CXX) $(CXXFLAGS) -O3 -g -DOMPI_SKIP_MPICXX -DNDEBUG $(SRC) -o $@

clean:
	rm -f step3_mpi