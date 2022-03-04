#-------------------------------------------------------------------------------
# KLU Lib/Makefile
#-------------------------------------------------------------------------------

LIBRARY = libklu
VERSION = 1.3.8
SO_VERSION = 1

default: static


include ../../SuiteSparse_config/SuiteSparse_config.mk

# compile and install in SuiteSparse/lib
library:
	$(MAKE) install INSTALL=$(SUITESPARSE)

# for testing only:
# TEST = -DTESTING

C = $(CC) $(CF)

INC = ../Include/klu_internal.h ../../SuiteSparse_config/SuiteSparse_config.h

I = -I../../AMD/Include -I../../BTF/Include -I../Include -I../../SuiteSparse_config


COMMON = klu_free_symbolic.o klu_defaults.o klu_analyze.o klu_memory.o

OBJ = $(COMMON)

static: $(AR_TARGET)

$(AR_TARGET): $(OBJ)
	@mkdir -p $(INSTALL_LIB)
	$(ARCHIVE) $(INSTALL_LIB)/$@ $^

$(OBJ): $(INC)

#-------------------------------------------------------------------------------

klu_analyze.o: ../klu_analyze.c
	$(C) -c $(I) $< -o $@

klu_defaults.o: ../klu_defaults.c
	$(C) -c $(I) $< -o $@

klu_free_symbolic.o: ../klu_free_symbolic.c
	$(C) -c $(I) $< -o $@

klu_memory.o: ../klu_memory.c
	$(C) -c $(I) $< -o $@

#-------------------------------------------------------------------------------

purge: distclean

distclean: clean
	- $(RM) -r $(PURGE)

clean:
	- $(RM) -r $(CLEAN)