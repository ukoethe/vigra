include $(VIGRA_PATH)/config/Makefile.include

CXXFLAGS = -O -I../include -I/software/fftw-2.1.3/include -I/home/rezniche/vigra/boost_test
LIB = -lm 
		 
test:: test.o
	$(CXX) $(CXXFLAGS) $(CXXINC) -o test test.o $(LIB)
	./test || { rm ./test; exit 1; }

clean::
	rm -f *.o test core res.xv *.plg *ilk *.d
	rm -rf Debug ii_files
	

%.d: %.cxx
	dependHeader.csh "$(CXXFLAGS)" "$(CXXINC)" $*.cxx 

%.h:
	@rm -f *.d

%.hxx:
	@rm -f *.d

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif
