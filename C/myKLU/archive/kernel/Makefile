CC=clang-13
CXX=clang++-13

CLIB=-L../../../KLU/lib -lklu -lm
INCL=-I../../../KLU/include

CHOLMOD_LIB= -lcholmod


STD=-std=c++17 -stdlib=libc++



klu_kernelx: klu_kernel.cpp Makefile
	$(CXX) $(INCL) $(STD) klu_kernel.cpp -o klu_kernel $(CLIB) -O3
	./klu_kernel

