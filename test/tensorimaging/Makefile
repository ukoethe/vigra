test: test.o
	$(CXX) -o test test.o
	./test || { rm ./test ; exit 1; }

tensorimaging.exe: test.d
	devenv "..\test.sln" /build Release /project tensorimaging
	./tensorimaging.exe || { rm ./tensorimaging.exe ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
