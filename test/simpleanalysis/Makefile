curdir=.

test: test.o
	$(CXX) -o test test.o
	./test || { rm ./test ; exit 1; }

simpleanalysis.exe: test.d
	devenv "..\test.sln" /build Release /project simpleanalysis
	./simpleanalysis.exe || { rm ./simpleanalysis.exe ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
