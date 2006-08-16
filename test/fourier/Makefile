curdir=.

TESTIMAGES = ghouse.gif filter.xv gaborresult.xv

all: test

include ../Makefile.include

test: test.o $(TESTIMAGES)
	$(CXX) -o test test.o $(IMPEX_LIB) $(FFTW_LIB)
	./test || { rm ./test$(EXEEXT) ; exit 1; }

fourier.exe: test.d
	devenv "..\test.sln" /build Release /project fourier
	./fourier.exe || { rm ./fourier.exe ; exit 1; }

$(TESTIMAGES):
	ln -s $(curdir)/$@ .

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

