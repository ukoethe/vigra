                 DOC = doxygen
	    DOCLATEX = doc++ -u -o $(DOCDIR)/latex/vigra.tex -t 
	      DOCDIR = ../doc
		 
FILES = \
	index.dxx 

all:: doc

doc::
	@if [ ! -d $(DOCDIR) ]; then mkdir $(DOCDIR); fi
	cat Doxyfile.in | sed "s/\(PROJECT_NUMBER *= \)/\1 $(VIGRA_VERSION)/" >! Doxyfile
	$(DOC)
	python post.py $(DOCDIR) $(VIGRA_VERSION)
	
latex::
	$(DOCLATEX)  $(FILES)
	
	
clean::
	rm -f $(DOCDIR)/*.html $(DOCDIR)/*.gif $(DOCDIR)/*.pnm \
              $(DOCDIR)/gifs.db $(DOCDIR)/*.ps $(DOCDIR)/*.eps \
              $(DOCDIR)/dxxgifs* $(DOCDIR)/*.class \
              $(DOCDIR)/formula.repository $(DOCDIR)/*.css
	
