test: test.o
	$(CXX) -o test test.o $(IMPEX_LIB)
	./test || { rm ./test ; exit 1; }

vc: image.exe
	./image.exe || { rm ./image.exe ; exit 1; }

image.exe:
	devenv "..\test.sln" /build Release /project image

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

EXTRA_CLEAN_TARGET=test.gif test.xv test.ras

include ../Makefile.include
