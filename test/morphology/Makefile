test: test.o
	$(CXX) -o test test.o
	./test || { rm ./test ; exit 1; }

vc: morphology.exe
	./morphology.exe || { rm ./morphology.exe ; exit 1; }

morphology.exe
	devenv "..\test.sln" /build Release /project morphology

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
