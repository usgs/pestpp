# This file is part of PEST++
top_builddir = ../..
include $(top_builddir)/global.mak

LIB := $(LIB_PRE)linear_analysis$(LIB_EXT)
OBJECTS := \
    covariance \
    linear_analysis \
    logger
OBJECTS := $(addsuffix $(OBJ_EXT),$(OBJECTS))


all: $(LIB)

$(LIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

clean:
	$(RM) $(OBJECTS) $(LIB)

.PHONY: all clean
