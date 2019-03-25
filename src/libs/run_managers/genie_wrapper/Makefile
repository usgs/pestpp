# This file is part of PEST++
top_builddir = ../../..
include $(top_builddir)/global.mak

LIB := $(LIB_PRE)rm_genie_wrapper$(LIB_EXT)
OBJECTS := RunManagerGenie$(OBJ_EXT)


all: $(LIB)

$(LIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

clean:
	$(RM) $(OBJECTS) $(LIB)

.PHONY: all clean
