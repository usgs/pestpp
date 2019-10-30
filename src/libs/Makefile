# This file is part of PEST++
top_builddir = ..
include $(top_builddir)/global.mak

SUBDIRS := \
    common \
    mio \
    pestpp_common \
    propack \
    run_managers \
    opt


all:	$(foreach d,$(SUBDIRS),$(d)-target)

clean:	$(foreach d,$(SUBDIRS),$(d)-clean)

.PHONY: all clean
