test: test.o
	$(CXX) -o test test.o
	./test || { rm ./test ; exit 1; }

error.exe: test.d
	devenv "..\test.sln" /build Release /project error
	./error.exe || { rm ./error.exe ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
