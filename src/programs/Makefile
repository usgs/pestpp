# This file is part of PEST++
top_builddir = ..
include $(top_builddir)/global.mak

SUBDIRS := \
    gsa \
    pest++ \
    pestpp-ies \
    pestpp-opt \
    sweep \
    pestpp-sqp \

# Does not seem to compile
# SUBDIRS += pestpp-ies

all:	$(foreach d,$(SUBDIRS),$(d)-target)

install:	$(foreach d,$(SUBDIRS),$(d)-install)

clean:	$(foreach d,$(SUBDIRS),$(d)-clean)

.PHONY: all install clean
