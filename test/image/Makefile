curdir=.

TESTIMAGES = lenna.xv

test: test.o $(TESTIMAGES)
	$(CXX) -o test test.o $(IMPEX_LIB)
	./test || { rm ./test ; exit 1; }

image.exe: test.d
	devenv "..\test.sln" /build Release /project image
	./image.exe || { rm ./image.exe ; exit 1; }

$(TESTIMAGES):
	ln -s $(curdir)/$@ .

EXTRA_CLEAN_TARGET=test.gif test.xv test.ras

include ../Makefile.include

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif
