include config/Makefile.include

all::
	@cd src ; $(MAKE) all ; cd ..

install::
	@cd src ; $(MAKE) install ; cd ..
	if test $(includedir) != "$(builddir)/include" ; then \
            $(INSTALL) -d $(includedir)/vigra ; \
            $(INSTALL) --mode=644 $(builddir)/include/vigra/*.hxx $(includedir)/vigra ; \
            $(INSTALL) --mode=644 $(builddir)/include/vigra/*.h $(includedir)/vigra ; \
    fi
	if test $(prefix) != $(builddir) ; then \
        $(INSTALL) -d $(prefix)/doc/documents ; \
            $(INSTALL) --mode=644 \
                 $(builddir)/doc/*.html \
                 $(builddir)/doc/classvigra*.gif $(builddir)/doc/form*.gif \
                 $(builddir)/doc/doxygen.gif $(builddir)/doc/doxygen.css \
               $(prefix)/doc ; \
            $(INSTALL) --mode=644 \
                 $(builddir)/doc/documents/*.ps \
                 $(builddir)/doc/documents/*.gif \
               $(prefix)/doc/documents ; \
    fi
	$(INSTALL) -d $(bindir)
	$(INSTALL) --mode=755 $(builddir)/vigra-config $(bindir)

examples::
	@cd src ; $(MAKE) examples ; cd ..

doc::
	cd docsrc; $(MAKE) VIGRA_VERSION=$(VIGRA_VERSION)

docclean::
	cd docsrc; $(MAKE) -i clean

clean::
	@cd src ; $(MAKE) clean ; cd ..

maintainer-clean:
	rm -f aclocal.m4 config.log config.cache config.status configure libtool
