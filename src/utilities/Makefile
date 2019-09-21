# This file is part of PEST++
top_builddir = ..
include $(top_builddir)/global.mak

SUBDIRS := \
    sweep


all:	$(foreach d,$(SUBDIRS),$(d)-target)

install:	$(foreach d,$(SUBDIRS),$(d)-install)

clean:	$(foreach d,$(SUBDIRS),$(d)-clean)

.PHONY: all install clean
