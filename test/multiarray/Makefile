curdir=.

test: test.o
	@test -d impex || mkdir impex
	$(CXX) -o test test.o $(IMPEX_LIB)
	./test || { rm ./test ; exit 1; }

multiarray.exe: test.d
	devenv "..\test.sln" /build Release /project multiarray
	./multiarray.exe || { rm ./multiarray.exe ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

EXTRA_CLEAN_TARGET=test?.png impex/test?.png

include ../Makefile.include
