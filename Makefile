ifneq "$(MAKECMDGOALS)" "autoconf"
include config/Makefile.include
endif

all::
	@cd src ; $(MAKE) all ; cd ..

install:: install-exec install-includes install-docs

install-exec:
	@cd src ; $(MAKE) install ; cd ..
	$(INSTALL) -d $(bindir)
	$(INSTALL) --mode=755 $(vigra_builddir)/vigra-config $(bindir)

install-includes:
	if test $(includedir) != "$(abs_top_builddir)/include" ; then \
            $(INSTALL) -d $(includedir)/vigra ; \
            $(INSTALL) --mode=644 $(abs_top_builddir)/include/vigra/*.hxx $(includedir)/vigra ; \
            $(INSTALL) --mode=644 $(abs_top_builddir)/include/vigra/*.h $(includedir)/vigra ; \
    fi

install-docs:
	if test $(prefix) != $(abs_top_builddir) ; then \
        $(INSTALL) -d $(prefix)/doc/documents ; \
            $(INSTALL) --mode=644 \
                 $(abs_top_builddir)/doc/*.html \
                 $(abs_top_builddir)/doc/classvigra*.gif $(abs_top_builddir)/doc/form*.gif \
                 $(abs_top_builddir)/doc/doxygen.gif $(abs_top_builddir)/doc/doxygen.css \
               $(prefix)/doc ; \
            $(INSTALL) --mode=644 \
                 $(abs_top_builddir)/doc/documents/*.ps \
                 $(abs_top_builddir)/doc/documents/*.gif \
               $(prefix)/doc/documents ; \
    fi

examples::
	@cd src ; $(MAKE) examples ; cd ..

doc::
	cd docsrc; $(MAKE) VIGRA_VERSION=$(VIGRA_VERSION)

docclean::
	cd docsrc; $(MAKE) -i clean

clean::
	@cd src ; $(MAKE) clean ; cd ..

distclean: clean
	rm -f config/vigra-config config/Makefile.include
	rm -f config.log config.cache config.status
	cp config/Makefile.include.empty config/Makefile.include

maintainer-clean: distclean
	rm -f aclocal.m4 configure libtool

autoconf:
	cd config && \
	aclocal --acdir=. && \
	autoconf && mv configure ..
#	autoconf --output=../configure
