test: test.o
	$(CXX) -o test test.o $(VIGRAIMPEX_LIB)
	./test || { rm ./test ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
