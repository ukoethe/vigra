test: test.o
	$(CXX) -o test test.o $(IMPEX_LIB)
	./test || { rm ./test ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

EXTRA_CLEAN_TARGET=test.gif test.xv test.ras

include ../Makefile.include
