# This file is part of PEST++
top_builddir = ../..
include $(top_builddir)/global.mak

EXE := pestpp-opt$(EXE_EXT)
OBJECTS := pestpp-opt$(OBJ_EXT)

PESTPP_INCLUDE += -I $(LIBS_DIR)/opt
PESTPP_LIBS := -L$(LIBS_DIR)/opt -lopt $(PESTPP_LIBS)


all: $(EXE)

$(EXE): $(OBJECTS)
	$(LD) $(LDFLAGS) $^ $(PESTPP_LIBS) -o $@

install: $(EXE)
	$(MKDIR) $(bindir)
	$(CP) $< $(bindir)

clean:
	$(RM) $(OBJECTS) $(EXE)

.PHONY: all install clean
