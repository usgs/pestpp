# This file is part of PEST++
top_builddir = ../..
include $(top_builddir)/global.mak

LIB := $(LIB_PRE)propack$(LIB_EXT)
OBJECTS := \
    dgemm_ovwr \
    dgetu0 \
    dlanbpro \
    dlansvd \
    dmgs.pentium \
    dreorth \
    dritzvec \
    dsafescal \
    propack_misc
OBJECTS := $(addsuffix $(OBJ_EXT),$(OBJECTS))


all: $(LIB)

$(LIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

clean:
	$(RM) $(OBJECTS) $(LIB)

.PHONY: all clean
