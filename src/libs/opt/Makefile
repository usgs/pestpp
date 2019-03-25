# This file is part of PEST++
top_builddir = ../..
include $(top_builddir)/global.mak

LIB := $(LIB_PRE)opt$(LIB_EXT)

# Compile all source files
SRC := $(wildcard *.cpp)
OBJECTS := $(SRC:.cpp=$(OBJ_EXT))

CPPFLAGS += -D HAVE_CMATH -D HAVE_FLOAT_H -D USE_EKKWSSMP


all: $(LIB)

# Build library on everything except ClpMain.cpp, which has main()
$(LIB): $(filter-out ClpMain$(OBJ_EXT),$(OBJECTS))
	$(AR) $(ARFLAGS) $@ $^

clean:
	$(RM) $(OBJECTS) $(LIB)

.PHONY: all clean
