all: test speedtest

test: test.o
	$(CXX) -o test test.o
	./test || { rm ./test ; exit 1; }

speedtest: speedtest.o
	$(CXX) -o speedtest speedtest.o
	./speedtest #|| { rm ./speedtest ; exit 1; }

multiconvolution.exe: test.d
	devenv "..\test.sln" /build Release /project multiconvolution
	./multiconvolution.exe || { rm ./multiconvolution.exe ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
