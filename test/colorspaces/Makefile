test: test.o
	$(CXX) $(CXXFLAGS) -o test test.o
	./test || { rm ./test ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
