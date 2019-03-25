# This file is part of PEST++
top_builddir = ../..
include $(top_builddir)/global.mak

EXE := pbin2ascii$(EXE_EXT)
OBJECTS := pbin2asciicpp$(OBJ_EXT)


all: $(EXE)

$(EXE): $(OBJECTS)
	$(LD) $(LDFLAGS) $^ $(PESTPP_LIBS) -o $@

install: $(EXE)
	$(MKDIR) $(bindir)
	$(CP) $< $(bindir)

clean:
	$(RM) $(OBJECTS) $(EXE)

.PHONY: all install clean
