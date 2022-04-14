CC=clang-13
CXX=clang++-13

CLIB=-L lib -lklu -lamd -lbtf -lsuitesparseconfig -lm
INCL=-I include

STD=-std=c++20 -stdlib=libc++

SUF=-Wall -O3

library:
	(cd SuiteSparse_config ; $(MAKE))
	(cd AMD/Lib ; $(MAKE))
	(cd BTF/Lib ; $(MAKE))
	(cd KLU/Lib ; $(MAKE))

klu_kernel: klu_kernel.cpp library
	$(CXX) $(INCL) $(STD) klu_kernel.cpp -o klu_kernel $(CLIB) $(SUF)
	./klu_kernel