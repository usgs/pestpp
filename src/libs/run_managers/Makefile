# This file is part of PEST++
top_builddir = ../..
include $(top_builddir)/global.mak

SUBDIRS := \
    abstract_base \
    external \
    serial \
    wrappers \
    yamr

all:	$(foreach d,$(SUBDIRS),$(d)-target)

clean:	$(foreach d,$(SUBDIRS),$(d)-clean)

.PHONY: all clean
