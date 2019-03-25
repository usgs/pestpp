# This file is part of PEST++
top_builddir = ../../..
include $(top_builddir)/global.mak

LIB := $(LIB_PRE)rm_genie$(LIB_EXT)
OBJECTS := \
    buffer \
    client \
    commandline \
    executable \
    GENIE_INTERFACE_BASE \
    genie_misc_routines \
    message \
    modelresult \
    modelrun \
    node \
    receivor \
    runprocess \
    socket_utilities \
    thread
OBJECTS := $(addsuffix $(OBJ_EXT),$(OBJECTS))


all: $(LIB)

$(LIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

clean:
	$(RM) $(OBJECTS) $(LIB)

.PHONY: all clean
