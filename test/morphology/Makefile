curdir=.

test: test.o
	$(CXX) -o test test.o
	./test || { rm ./test ; exit 1; }

morphology.exe: test.d
	devenv "..\test.sln" /build Release /project morphology
	./morphology.exe || { rm ./morphology.exe ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
