test: test.o
	$(CXX) -o test test.o
	./test || { rm ./test ; exit 1; }

vc: simpleanalysis.exe
	./simpleanalysis.exe || { rm ./simpleanalysis.exe ; exit 1; }

simpleanalysis.exe
	devenv "..\test.sln" /build Release /project simpleanalysis

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
