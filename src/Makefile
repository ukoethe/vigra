SUBDIRS = \
        impex \
        examples

all::
	cd impex; $(MAKE)
	cd examples; $(MAKE)

clean::
	rm -f *.ncb *.opt 
	cd impex; $(MAKE) -i clean
	cd examples; $(MAKE) -i clean
      
