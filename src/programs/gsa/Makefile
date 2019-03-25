# This file is part of PEST++
top_builddir = ../..
include $(top_builddir)/global.mak

EXE := pestpp-sen$(EXE_EXT)
OBJECTS := \
    main \
    GsaAbstractBase \
    MorrisMethod \
    PooledVariance \
    sobol \
    Stats
OBJECTS := $(addsuffix $(OBJ_EXT),$(OBJECTS))


all: $(EXE)

$(EXE): $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $^ $(PESTPP_LIBS)

install: $(EXE)
	$(MKDIR) $(bindir)
	$(CP) $< $(bindir)

clean:
	$(RM) $(OBJECTS) $(EXE)

.PHONY: all install clean
