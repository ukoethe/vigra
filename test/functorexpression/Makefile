include $(VIGRA_PATH)/config/Makefile.include

		   CXXFLAGS = -O -I../include
	                LIB = -lm
		 
test:: test.o
	$(CXX) $(CXXFLAGS) $(CXXINC) -o test test.o $(LIB)
	./test || { rm ./test; exit 1; }

clean::
	rm -f *.o test core res.xv *.plg *ilk *.d
	rm -rf Debug ii_files
	

%.d: %.cxx
	dependHeader.csh "$(CXXFLAGS)" "$(CXXINC)" $*.cxx 

include test.d
