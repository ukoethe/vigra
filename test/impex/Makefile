include $(VIGRA_PATH)/config/Makefile.include

CXXFLAGS = -O -I../include # -I../../experiments/newimpex/include
LIB = -lvigraimpex -ljpeg -ltiff -lm

CXXLIBDIR = \
 -L/software/libpng-1.2.1/SunOS-5.8/lib \
 -L/software/jpeg-6b/SunOS-5.6/lib \
 -L/software/tiff-v3.5.5/SunOS-5.6/lib \
 -L../../current/lib/Solaris
 
LIB = -lpng -ltiff -ljpeg -lvigraimpex
                 		   
test:  test.o ../../current/lib/Solaris/libvigraimpex.a
	$(CXX) $(CXXFLAGS) $(CXXINC) -o test test.o $(CXXLIBDIR) $(LIB)
#	./test || { rm ./test ; exit 1; }

clean:
	rm -f *.o test core res.* *.plg *ilk *.d
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
