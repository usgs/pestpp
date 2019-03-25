# This file is part of PEST++
top_builddir = ../..
include $(top_builddir)/global.mak

LIB := $(LIB_PRE)mio$(LIB_EXT)
OBJECTS := mio$(OBJ_EXT)

FFLAGS += $(FFREE)


all: $(LIB)

$(LIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

clean:
	$(RM) $(OBJECTS) $(LIB) *.mod

.PHONY: all clean
