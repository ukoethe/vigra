curdir=.

test: test.o
#	$(CXX) -o test test.o $(IMPEX_LIB) $(FFTW3_LIB)
	$(CXX) -o test test.o $(IMPEX_LIB) -L/software/fftw-3.0.1/SunOS-5.8/lib -lfftw3 -lm -Wl,--rpath -Wl,/software/fftw-3.0.1/SunOS-5.8/lib
	./test || { rm ./test ; exit 1; }

fourier.exe: test.d
	devenv "..\test.sln" /build Release /project fourier
	./fourier.exe || { rm ./fourier.exe ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
CPPFLAGS = -I/home/koethe/vigra/current/include -I../include -I/software/fftw-3.0.1/include
