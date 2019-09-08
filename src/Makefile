# Top-level GNU Make file for PEST++

# Customize options in global.mak
top_builddir = .
include $(top_builddir)/global.mak


all: libs-target programs-target utilities-target

libs-target:
	(cd libs; $(MAKE))
	touch libs/build-stamp

programs-target: libs-target
	(cd programs; $(MAKE))

install: libs-target
	(cd programs; $(MAKE) install)

clean:
	(cd programs; $(MAKE) clean)
	(cd libs; $(MAKE) clean)
	$(RM) libs/build-stamp

.PHONY: all libs-target programs-target utilities-target install clean
