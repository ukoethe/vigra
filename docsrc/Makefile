                 DOC = $(HOME)/C++/doc++/src/doc++ -G -j -u -d $(DOCDIR) -B footer.html
	    DOCLATEX = doc++ -u -o $(DOCDIR)/latex/vigra.tex -t 
	      DOCDIR = ../doc
		 
FILES = \
	index.dxx 

all:: doc

doc::
	@if [ ! -d $(DOCDIR) ]; then mkdir $(DOCDIR); fi
	@makefooter.csh $(VIGRA_VERSION)
	$(DOC) $(FILES)
	python postprocessing.py $(DOCDIR) $(VIGRA_VERSION)
	
latex::
	$(DOCLATEX)  $(FILES)
	
	
clean::
	rm -f $(DOCDIR)/*.html $(DOCDIR)/*.gif $(DOCDIR)/*.pnm \
              $(DOCDIR)/gifs.db $(DOCDIR)/*.ps $(DOCDIR)/*.eps \
              $(DOCDIR)/dxxgifs* $(DOCDIR)/*.class
	
