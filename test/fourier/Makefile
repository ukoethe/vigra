test: test.o
#	$(CXX) -o test test.o $(IMPEX_LIB) $(FFTW3_LIB)
	$(CXX) -o test test.o $(IMPEX_LIB) /software/fftw-3.0.1/fftw-3.0.1/.libs/libfftw3.a -lm
	./test || { rm ./test ; exit 1; }

fourier.exe: test.d
	devenv "..\test.sln" /build Release /project fourier
	./fourier.exe || { rm ./fourier.exe ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
CPPFLAGS = -I/home/koethe/vigra/current/include -I../include -I/software/fftw-3.0.1/fftw-3.0.1/api
