test: test.o
	$(CXX) -o test test.o
	./test || { rm ./test ; exit 1; }

vc: error.exe
	./error.exe || { rm ./error.exe ; exit 1; }

error.exe:
	devenv "..\test.sln" /build Release /project error

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
