CC=clang-13
CXX=clang++-13

CLIB=-L../KLU/lib -lklu -lbtf -lamd -lcolamd -lsuitesparseconfig -lccolamd -lcamd -lmetis -lm -lrt -lcholmod -llapack -lblas

INCL=-I../KLU/include

STD=-std=c++17 -stdlib=libc++

sparse3: sparse3.cpp Makefile
	$(CXX) $(INCL) $(STD) sparse3.cpp -o sparse3 $(CLIB) -O2
	- ./sparse3 < ../Matrix_Sample/test2.mtx