curdir=.

test: test.o
	$(CXX) -o test test.o $(IMPEX_LIB)
	./test || { rm ./test ; exit 1; }

image.exe: test.d
	devenv "..\test.sln" /build Release /project utilities
	./utilities.exe || { rm ./utilities.exe ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

EXTRA_CLEAN_TARGET=

include ../Makefile.include
