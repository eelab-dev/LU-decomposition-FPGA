#-------------------------------------------------------------------------------
# SuiteSparse_config Makefile
#-------------------------------------------------------------------------------

SUITESPARSE ?= $(realpath $(CURDIR)/..)
export SUITESPARSE

# version of SuiteSparse_config is also version of SuiteSparse meta-package
LIBRARY = libsuitesparseconfig
VERSION = 5.10.1
SO_VERSION = 5

default: static

include SuiteSparse_config.mk


OBJ = SuiteSparse_config.o

SuiteSparse_config.o: SuiteSparse_config.c SuiteSparse_config.h
	$(CC) $(CF) -c SuiteSparse_config.c

static: $(AR_TARGET)

$(AR_TARGET): $(OBJ)
	$(ARCHIVE) $(AR_TARGET) SuiteSparse_config.o
	$(CP) $@ $(INSTALL_LIB)/$@

distclean: purge

purge: clean
	( cd xerbla ; $(MAKE) purge )
	- $(RM) -r $(PURGE)

clean:
	( cd xerbla ; $(MAKE) clean )
	- $(RM) -r $(CLEAN)

