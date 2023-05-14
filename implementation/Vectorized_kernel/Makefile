## This file is used to compile the vectorized code

CXX ?= g++
CXXFLAGS = -O3 -Wall -Wextra -pedantic -std=c++11

.PHONY: clean

main: main.cpp mp2_reduction_kernel.o mp2_reduction_kernel_avx.o
	$(CXX) $(CXXFLAGS) -lopenblas -D_ISPC_ -o $@ $^

mp2_reduction_kernel.o: mp2_reduction_kernel.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
mp2_reduction_kernel_avx.o: mp2_reduction_kernel.ispc
	ispc -O3 --arch=x86-64 --target=avx2-i64x4 -o $@ $<
clean:
	rm -f main mp2_reduction_kernel*.o