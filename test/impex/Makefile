include $(VIGRA_PATH)/config/Makefile.include

                   CXXFLAGS = -O -I../include 
                        LIB = -lvigraimpex -ljpeg -ltiff -lm
                 		   
test:  test.o 
	$(CXX) $(CXXFLAGS) $(CXXINC) -o test test.o $(CXXLIBDIR) $(LIB)
	./test || { rm ./test ; exit 1; }

clean:
	rm -f *.o test core res.* *.plg *ilk *.d
	rm -rf Debug ii_files

%.d: %.cxx
	dependHeader.csh "$(CXXFLAGS)" "$(CXXINC)" $*.cxx 

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif
