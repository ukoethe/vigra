curdir=.

TESTIMAGES = lenna.xv  lenna_gifref.xv  lennafloat.xv  lennafloatrgb.xv  lennargb.xv

all: test

include ../Makefile.include

CXXDEBUG := -O # or -g

test: test.o $(vigra_builddir)/src/impex/.libs/libvigraimpex.a $(TESTIMAGES)
	$(CXX) -o test test.o $(vigra_builddir)/src/impex/.libs/libvigraimpex.a $(IMPEX_DEPLIBS)
	./test $(curdir) || { rm ./test ; exit 1; }

impex.exe: test.d
	devenv "..\test.sln" /build Release /project impex
	./impex.exe || { rm ./impex.exe ; exit 1; }

$(TESTIMAGES):
	ln -s $(curdir)/$@ .

EXTRA_CLEAN_TARGET = res.* resrgb.*

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

