test: test.o
	$(CXX) -o test test.o $(IMPEX_LIB)
	./test || { rm ./test ; exit 1; }

convolution.exe: test.d
	devenv "..\test.sln" /build Release /project convolution
	./convolution.exe || { rm ./convolution.exe ; exit 1; }
        
ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
