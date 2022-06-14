CC=clang
CXX=clang++

CLIB=-L../../SuiteSparse/lib -lklu -lbtf -lamd -lcolamd -lsuitesparseconfig -lccolamd -lcamd -lmetis -lm -lrt -lcholmod -llapack -lblas

CXLIB=-L../../SuiteSparse/lib -lklu

INCL=-I../../SuiteSparse/include

STD=-std=c++17 -stdlib=libc++

LIBS=-L../../SuiteSparse/lib -lumfpack

klu_simple: klu_simple.cpp Makefile
	$(CXX) $(INCL) $(STD) klu_simple.cpp -o klu_simple $(CXLIB)
	./klu_simple
