CC=clang-13
CXX=clang++-13

CLIB=-L../../KLU/lib -lklu -lm
INCL=-I../../KLU/include

CHOLMOD_LIB= -lcholmod


STD=-std=c++17 -stdlib=libc++

klu_kernel: klu_kernel.c Makefile
	$(CC) $(INCL) klu_kernel.c -o klu_kernel $(CLIB)
	# ./klu_kernel

klu_kernel2: klu_kernel2.c Makefile
	$(CC) $(INCL) klu_kernel2.c -o klu_kernel2 $(CLIB) $(CHOLMOD_LIB)
	# ./klu_kernel2 <

klu_kernel3: klu_kernel3.c Makefile
	$(CC) $(INCL) klu_kernel3.c -o klu_kernel3 $(CLIB) $(CHOLMOD_LIB)
	# ./klu_kernel2 <

klu_kernel4: klu_kernel4.c Makefile
	$(CC) $(INCL) klu_kernel4.c -o klu_kernel4 $(CLIB)
	# ./klu_kernel4 <

klu_kernel5: klu_kernel5.cpp Makefile
	$(CXX) $(INCL) $(STD) klu_kernel5.cpp -o klu_kernel5 $(CLIB) -O3
	./klu_kernel5

klu_kernelx: klu_kernel.cpp Makefile
	$(CXX) $(INCL) $(STD) klu_kernel.cpp -o klu_kernel $(CLIB) -O3
	./klu_kernel

klu_kernelx2: klu_kernel2.cpp Makefile
	$(CXX) $(INCL) $(STD) klu_kernel2.cpp -o klu_kernel $(CLIB) -O3
	./klu_kernel