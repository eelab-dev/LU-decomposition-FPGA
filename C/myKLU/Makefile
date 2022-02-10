CC=clang-13
CXX=clang++-13

CLIB=-L../../KLU/lib -lklu -lm
INCL=-I../../KLU/include

CHOLMOD_LIB= -lcholmod


STD=-std=c++17 -stdlib=libc++

klu_kernel: klu_kernel.c Makefile
	$(CC) $(INCL) klu_kernel.c -o klu_kernel -lm
	# ./klu_kernel

klu_kernel2: klu_kernel2.c Makefile
	$(CC) $(INCL) klu_kernel2.c -o klu_kernel2 $(CLIB) $(CHOLMOD_LIB)
	# ./klu_kernel2 <