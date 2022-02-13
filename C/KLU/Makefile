CC=clang-13
CXX=clang++-13

CLIB=-L../../KLU/lib -lklu -lbtf -lamd -lcolamd -lsuitesparseconfig -lccolamd -lcamd -lmetis -lm -lrt -lcholmod -llapack -lblas

INCL=-I../../KLU/include

STD=-std=c++17 -stdlib=libc++

klu_simple: klu_simple.c Makefile
	$(CC) $(INCL) klu_simple.c -o klu_simple $(CLIB)
	./klu_simple

klu_demo: kludemo.c Makefile
	$(CC) $(INCL) kludemo.c -o kludemo $(CLIB)
	- ./kludemo < ../../Matrix_Sample/1c.mtx
	- ./kludemo < ../../Matrix_Sample/arrowc.mtx
	- ./kludemo < ../../Matrix_Sample/arrow.mtx
	- ./kludemo < ../../Matrix_Sample/impcol_a.mtx
	- ./kludemo < ../../Matrix_Sample/w156.mtx
	- ./kludemo < ../../Matrix_Sample/ctina.mtx

klu_demo2: kluldemo.cpp Makefile
	$(CXX) $(INCL) $(STD) kluldemo.cpp -o kluldemo $(CLIB)
	- ./kluldemo < ../../Matrix_Sample/bcsstk06.mtx

# klu_extract: klu_extract.cpp Makefile
# 	$(CXX) $(INCL) $(STD) klu_extract.cpp -o klu_extract $(CLIB)
# 	./klu_extract < ../../Matrix_Sample/test.mtx

klu_extract: klu_extract.c Makefile
	$(CC) $(INCL) klu_extract.c -o klu_extract $(CLIB)
	./klu_extract < ../../Matrix_Sample/test.mtx

klu_solve: klu_solve.c Makefile
	$(CC) $(INCL) klu_solve.c -o klu_solve $(CLIB)

btf: btf.c Makefile
	$(CC) $(INCL) btf.c -o btf $(CLIB)
	./btf < ../../Matrix_Sample/test.mtx