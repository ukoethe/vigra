SUBDIRS = \
        src 

all::
	cd src; $(MAKE)

doc::
	cd docsrc; $(MAKE) VIGRA_VERSION=$(VIGRA_VERSION)

clean::
	cd src; $(MAKE) -i clean
      
docclean::
	cd docsrc; $(MAKE) -i clean

realclean:: clean
	rm -rf lib/*

	
