#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8

SUBDIRS := ../../libana .

lib_LTLIBRARIES = libpython_interface.la

libpython_interface_la_SOURCES = \
        your_module.fpp

EXTRA_DIST = \
	Makefile.depends 

libpython_interface_la_LDFLAGS = -module -avoid-version -export-dynamic
libpython_interface_la_LIBADD = -llapack -lblas ../../../lib/lib.a ../../libana/libana.a
libpython_interface_la_FCFLAGS = $(FCFLAGS) $(PYTHON_CPPFLAGS)

# Add these lines
CCLD = $(FC)
LINK = $(LIBTOOL) --tag=FC --mode=link $(CCLD) $(AM_CFLAGS) $(CFLAGS) \
       $(AM_LDFLAGS) $(LDFLAGS) -o $@

.fpp.lo:
if USEKCOMP
	cp $< $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(LIBTOOL) --tag=FC --mode=compile $(FC) $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana $(libpython_interface_la_FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $< $*.f90
	$(LIBTOOL) --tag=FC --mode=compile $(FC) $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana $(libpython_interface_la_FCFLAGS) -c $*.f90
endif

# コンパイルフラグにPythonヘッダーのパスを追加
#AM_CPPFLAGS = $(PYTHON_CPPFLAGS)
#AM_LDFLAGS = -module -avoid-version -export-dynamic

#1. FCコンパイラを使用してFortranコードをコンパイルすることで、適切なフラグが使用
#2. リンク時にCCFLAGSではなくFCFLAGSとLDFLAGSを使用することで、不適切なフラグ（特に -g）が渡されるのを防いだ
#3. libtoolに --tag=FC オプションを使用することで、Fortranコードの適切な処理が保証された

#libpython_interface_la_LDFLAGS = -module -avoid-version -export-dynamic -O2 -fopenmp
#libpython_interface_la_CPPFLAGS = $(PYTHON_CPPFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
