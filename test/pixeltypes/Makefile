test: test.o
	$(CXX) -o test test.o
	./test || { rm ./test ; exit 1; }

vc: pixeltypes.exe
	./pixeltypes.exe || { rm ./pixeltypes.exe ; exit 1; }

pixeltypes.exe
	devenv "..\test.sln" /build Release /project pixeltypes

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
