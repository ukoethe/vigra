include $(VIGRA_PATH)/config/Makefile.include

		   CXXFLAGS = -O -I../include 
	                LIB = -lvigraimpex -ltiff -ljpeg
                 
test:  test.o 
	$(CXX) $(CXXFLAGS) $(CXXINC) -o test test.o $(CXXLIBDIR) $(LIB)
	./test || { rm ./test ; exit 1; }
        
clean::
	rm -f *.o test core test.gif test.xv test.ras *.plg *ilk *.d
	rm -rf Debug ii_files
	
	
%.d: %.cxx
	dependHeader.csh "$(CXXFLAGS)" "$(CXXINC)" $*.cxx 

include test.d
