#-------------------------------------------------------------------------------
# BTF Lib/Makefile
#-------------------------------------------------------------------------------

LIBRARY = libbtf
VERSION = 1.2.6
SO_VERSION = 1

default: static

include ../../SuiteSparse_config/SuiteSparse_config.mk

# for testing only:
# TEST = -DTESTING

C = $(CC) $(CF)

INC = ../Include/btf.h ../Include/btf_internal.h

I = -I../Include -I../../SuiteSparse_config

all: static

OBJ = btf_order.o btf_maxtrans.o btf_strongcomp.o

static: $(AR_TARGET)

$(AR_TARGET): $(OBJ)
	$(ARCHIVE) $@ $^
	$(CP) $@ $(INSTALL_LIB)/$@

#-------------------------------------------------------------------------------

btf_order.o: ../btf_order.c
	$(C) -c $(I) $< -o $@

btf_maxtrans.o: ../btf_maxtrans.c
	$(C) -c $(I) $< -o $@

btf_strongcomp.o: ../btf_strongcomp.c
	$(C) -c $(I) $< -o $@

#-------------------------------------------------------------------------------

purge: distclean

distclean: clean
	- $(RM) -r $(PURGE)

clean:
	- $(RM) -r $(CLEAN)
