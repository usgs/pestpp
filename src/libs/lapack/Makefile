#
#  Top Level Makefile for LAPACK
#  Version 3.4.1
#  April 2012
#

include make.inc

all: lib blaslib

lib: lapacklib
#lib: blaslib variants lapacklib tmglib

lapack_install:
	$(MAKE) -C INSTALL run

blaslib:
	$(MAKE) -C BLAS

lapacklib:
	$(MAKE) -C SRC double

clean:
	$(MAKE) -C INSTALL clean
	$(MAKE) -C BLAS clean
	$(MAKE) -C SRC clean
	rm -f *.a
