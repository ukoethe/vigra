test: test.o
	$(CXX) -o test test.o
	./test || { rm ./test ; exit 1; }

vc: pixeltypes.exe
	devenv "..\test.sln" /build Release /project pixeltypes
	./pixeltypes.exe || { rm ./pixeltypes.exe ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
