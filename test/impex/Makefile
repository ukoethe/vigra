test: test.o
	$(CXX) -o test test.o $(vigra_builddir)/src/impex/libvigraimpex.a $(IMPEX_DEPLIBS)
	./test || { rm ./test ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
