                 DOC = doc++ -G -j -u -d $(DOCDIR) -B footer.html
	    DOCLATEX = doc++ -u -o $(DOCDIR)/latex/vigra.tex -t 
	      DOCDIR = ../doc
		 
FILES = \
	index.dxx 

all:: doc

doc::
	@if [ ! -d $(DOCDIR) ]; then mkdir $(DOCDIR); fi
	@makefooter.csh $(VIGRA_VERSION)
	$(DOC) $(FILES)
	@for i in $(DOCDIR)/*.html ; do \
		echo $$i; \
		cat $$i | sed -f postprocess_doc.sed > /tmp/t.html ; \
		mv /tmp/t.html $$i; \
	done
	@cat $(DOCDIR)/aindex.html | \
              sed 's?<H2>Variables</H2>?<H2>typedefs</H2>?' > /tmp/t.html; \
              mv /tmp/t.html $(DOCDIR)/aindex.html
	@cat $(DOCDIR)/index.html | \
              sed 's?VERSION_VERSION_VERSION?$(VIGRA_VERSION)?g' > \
              /tmp/t.html; mv /tmp/t.html $(DOCDIR)/index.html
	
latex::
	$(DOCLATEX)  $(FILES)
	
	
clean::
	rm -f $(DOCDIR)/*.html $(DOCDIR)/*.gif $(DOCDIR)/*.pnm \
              $(DOCDIR)/gifs.db $(DOCDIR)/*.ps $(DOCDIR)/*.eps \
              $(DOCDIR)/dxxgifs* $(DOCDIR)/*.class
	
