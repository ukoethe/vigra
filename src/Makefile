include ../config/Makefile.include

all::
	@cd impex ; $(MAKE) all ; cd ..

install::
	@cd impex ; $(MAKE) install ; cd ..

examples::
	@cd examples ; $(MAKE) CXX=$(CXX) all ; cd ..

clean::
	@cd impex ; $(MAKE) clean ; cd ..
	@cd examples ; $(MAKE) clean ; cd ..
