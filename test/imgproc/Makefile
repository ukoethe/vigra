test: test.o
	$(CXX) -o test test.o $(IMPEX_LIB)
	./test || { rm ./test ; exit 1; }

imgproc.exe: test.d
	devenv "..\test.sln" /build Release /project imgproc
	./imgproc.exe || { rm ./imgproc.exe ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
