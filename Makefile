ifneq "$(MAKECMDGOALS)" "autoconf"
include config/Makefile.include
endif

all::
	@cd src ; $(MAKE) all ; cd ..

install:: install-exec install-includes install-docs

install-exec:
	@cd src ; $(MAKE) install ; cd ..
	$(INSTALL) -d $(bindir)
	$(INSTALL) --mode=755 $(vigra_builddir)/config/vigra-config $(bindir)

install-includes:
	if test "$(includedir)" != "$(vigra_srcdir)/include" ; then \
          $(INSTALL) -d $(includedir)/vigra ; \
          $(INSTALL) --mode=644 $(vigra_srcdir)/include/vigra/*.hxx $(includedir)/vigra ; \
        fi
install-headers: install-includes

install-docs:
	$(INSTALL) -d $(docdir)
	$(INSTALL) --mode=644 LICENSE $(docdir)
	if test "$(docdir)" != "$(vigra_srcdir)/doc/vigra" ; then \
          $(INSTALL) --mode=644 \
            $(vigra_srcdir)/doc/vigra/*.html \
            $(vigra_srcdir)/doc/vigra/classvigra*.png $(vigra_srcdir)/doc/vigra/form*.png \
            $(vigra_srcdir)/doc/vigra/doxygen.png $(vigra_srcdir)/doc/vigra/doxygen.css \
            $(docdir) ; \
          $(INSTALL) -d $(docdir)/documents ; \
          $(INSTALL) --mode=644 \
            $(vigra_srcdir)/doc/vigra/documents/*.ps \
            $(vigra_srcdir)/doc/vigra/documents/*.gif \
            $(docdir)/documents ; \
        fi

examples::
	@cd src ; $(MAKE) examples ; cd ..

doc::
	cd docsrc; $(MAKE) VIGRA_VERSION=$(VIGRA_VERSION)

docclean::
	cd docsrc; $(MAKE) -i clean

clean::
	@cd src ; $(MAKE) clean ; cd ..
	rm -f a.exe a.out

distclean: clean
	rm -f config/vigra-config config/Makefile.include
	rm -f config.log config.cache config.status libtool
	cp config/Makefile.include.empty config/Makefile.include

maintainer-clean: distclean docclean
	rm -f config/aclocal.m4 configure

autoconf:
	cd config && \
	aclocal --acdir=. && \
	autoconf && mv configure ..

libtoolize:
	@echo "*** SET $LTPARAMS TO EXTRA PARAMS (LIKE --force) IF YOU WISH ***"
	@libtool --version
	@test -h configure.in || ln -s config/configure.in
	@libtoolize -c $(LTPARAMS)
	@test -h configure.in && rm configure.in

.PHONY: autoconf all clean distclean doc docclean examples install install-docs install-exec install-includes install-headers
