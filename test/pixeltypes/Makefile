include $(VIGRA_PATH)/config/Makefile.include


CXXFLAGS =  -O -I../include
LIB = 
		 
test:: test.o
	$(CXX) $(CXXFLAGS) $(CXXINC) -o test test.o 
	./test || { rm ./test; exit 1; }
	
clean::
	rm -f *.o *.d *.bak graphtest

%.d: %.cxx
	dependHeader.csh "$(CXXFLAGS)" "$(CXXINC)" $*.cxx

%.h:
	@rm -f *.d

%.hxx:
	@rm -f *.d
	
ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif
