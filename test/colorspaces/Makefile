test: test.o
	$(CXX) $(CXXFLAGS) -o test test.o
	./test || { rm ./test ; exit 1; }

vc: colorspaces.exe
	./colorspaces.exe || { rm ./colorspaces.exe ; exit 1; }

colorspaces.exe:
	devenv "..\test.sln" /build Release /project colorspaces

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
