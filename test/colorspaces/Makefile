test: test.o
	$(CXX) $(CXXFLAGS) -o test test.o
	./test || { rm ./test ; exit 1; }

colorspaces.exe: test.d
	devenv "..\test.sln" /build Release /project colorspaces
	./colorspaces.exe || { rm ./colorspaces.exe ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
