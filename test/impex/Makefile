curdir=.
ifneq "$(MAKECMDGOALS)" "clean"
include ../../current/build.cygwin.gcc3.4/config/Makefile.include
endif

CPPFLAGS += -I../include

#EXTRACXXFLAGS =  -DNEWIMPEX
#EXTRACXXINCLUDES =  -I../../experiments/newimpex/include
CXXDEBUG = -g # or -O

#LIBS = -L../../experiments/newimpex/src/impex  -L/software/libpng-1.2.1/SunOS-5.8/lib -lpng $(VIGRAIMPEX_LIB)

test: test.o $(vigra_builddir)/src/impex/.libs/libvigraimpex.a
	$(CXX) -o test test.o $(vigra_builddir)/src/impex/.libs/libvigraimpex.a $(IMPEX_DEPLIBS)
	./test || { rm ./test ; exit 1; }

impex.exe: test.d
	devenv "..\test.sln" /build Release /project impex
	./impex.exe || { rm ./impex.exe ; exit 1; }

clean::
	rm -f *.o test core *.d $(EXTRA_CLEAN_TARGET)
	rm -rf Debug ii_files *.plg *.ilk

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

