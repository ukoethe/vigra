SUBDIRS = \
        impex \
        examples

all:: libs examples

libs::
	cd impex; $(MAKE)

examples::
	cd examples; $(MAKE)

clean::
	rm -f *.ncb *.opt 
	cd impex; $(MAKE) -i clean
	cd examples; $(MAKE) -i clean
