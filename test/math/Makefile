test: test.o
	$(CXX) -o test test.o
	./test || { rm ./test ; exit 1; }

math.exe: test.d
	devenv "..\test.sln" /build Release /project math
	./math.exe || { rm ./math.exe ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
