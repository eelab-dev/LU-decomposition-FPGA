CC=clang
CXX=clang++

CLIB=-L lib -lklu -lamd -lbtf -lsuitesparseconfig -lm
INCL=-I include -I myKLU

# STD=-std=c++20 -stdlib=libc++

SUF=-Wall -O3

library:
	(cd SuiteSparse_config ; $(MAKE))
	(cd AMD/Lib ; $(MAKE))
	(cd BTF/Lib ; $(MAKE))
	(cd KLU/Lib ; $(MAKE))

klu_kernel: klu_kernel.cpp library
	$(CXX) $(INCL) klu_kernel.cpp -o klu_kernel $(CLIB) $(SUF)
	./klu_kernel

klu_bench: klu_bench.cpp library
	$(CXX) $(INCL) klu_bench.cpp -o klu_bench $(CLIB) $(SUF)
	./klu_bench