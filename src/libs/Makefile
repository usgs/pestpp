# This file is part of PEST++
top_builddir = ..
include $(top_builddir)/global.mak

SUBDIRS := \
    common \
    linear_analysis \
    mio \
    pestpp_common \
    pest_routines \
    propack \
    run_managers \
    opt


all:	$(foreach d,$(SUBDIRS),$(d)-target)

clean:	$(foreach d,$(SUBDIRS),$(d)-clean)

.PHONY: all clean
