SUBDIRS = \
        src

all:: libs examples

libs::
	cd src; $(MAKE) libs

examples::
	cd src; $(MAKE) examples

doc::
	cd docsrc; $(MAKE) VIGRA_VERSION=$(VIGRA_VERSION)

clean::
	cd src; $(MAKE) -i clean

docclean::
	cd docsrc; $(MAKE) -i clean

realclean:: clean
	rm -rf lib/*
