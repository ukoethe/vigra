CXXFLAGS = -O -I../include `vigra-config --cxxflags`

test: test.o
	$(CXX) $(CXXFLAGS) -o test test.o `vigra-config --libs`
	./test || { rm ./test ; exit 1; }

clean:
	rm -f test.gif test.xv test.ras
	rm -f *.o test core *.d
	rm -rf Debug ii_files *.plg *.ilk

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $*.cxx

%.d: %.cxx
	@echo updating header dependencies for $*.cxx
	@rm -f $*.d
	@(echo "$*.d \c" ; $(CXX) -M $(CXXFLAGS) -c $*.cxx) > $*.d \
        || { rm -f $*.d; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif
