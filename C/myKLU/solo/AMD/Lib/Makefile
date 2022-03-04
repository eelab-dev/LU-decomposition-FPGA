#-------------------------------------------------------------------------------
# AMD Lib/Makefile
#-------------------------------------------------------------------------------

LIBRARY = libamd
VERSION = 2.4.6
SO_VERSION = 2

default: static

include ../../SuiteSparse_config/SuiteSparse_config.mk

C = $(CC) $(CF) -I../Include -I../../SuiteSparse_config

#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------

AMD = amd_aat amd_1 amd_2 amd_dump amd_postorder amd_defaults \
        amd_post_tree  \
	amd_order amd_control amd_info amd_valid amd_preprocess

INC = ../Include/amd.h ../Include/amd_internal.h \
      ../../SuiteSparse_config/SuiteSparse_config.h

#-------------------------------------------------------------------------------
# object files for each version
#-------------------------------------------------------------------------------

AMDI = $(addsuffix .o, $(subst amd_,amd_i_,$(AMD)))
AMDL = $(addsuffix .o, $(subst amd_,amd_l_,$(AMD)))
OBJ = $(AMDI) $(AMDL)

#-------------------------------------------------------------------------------
# compile each int and long routine (with no real/complex version)
#-------------------------------------------------------------------------------

amd_i_%.o: ../amd_%.c $(INC)
	$(C) -DDINT -c $< -o $@

amd_l_%.o: ../amd_%.c $(INC)
	$(C) -DDLONG -c $< -o $@

#-------------------------------------------------------------------------------
# Create the static library (C versions only)
#-------------------------------------------------------------------------------

static: $(AR_TARGET)

$(AR_TARGET): $(OBJ)
	@mkdir -p $(INSTALL_LIB)
	$(ARCHIVE) $(INSTALL_LIB)/$@ $^

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

clean:
	- $(RM) -r $(CLEAN)

purge: distclean

distclean: clean
	- $(RM) -r $(PURGE)

