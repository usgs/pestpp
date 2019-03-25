# This file is part of PEST++
top_builddir = ..
include $(top_builddir)/global.mak

SUBDIRS := \
    ascii2pbin \
    pbin2ascii \
    pbin_dump \
    sweep


all:	$(foreach d,$(SUBDIRS),$(d)-target)

install:	$(foreach d,$(SUBDIRS),$(d)-install)

clean:	$(foreach d,$(SUBDIRS),$(d)-clean)

.PHONY: all install clean
