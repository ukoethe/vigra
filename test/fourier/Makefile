curdir=.

all: test

include ../Makefile.include

test: test.o
	$(CXX) -o test test.o $(IMPEX_LIB) $(FFTW_LIB)
	./test || { rm ./test ; exit 1; }

fourier.exe: test.d
	devenv "..\test.sln" /build Release /project fourier
	./fourier.exe || { rm ./fourier.exe ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

