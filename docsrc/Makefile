         DOC = doxygen
    DOCLATEX = doc++ -u -o $(DOCDIR)/latex/vigra.tex -t
      DOCDIR = ../doc/vigra

FILES = \
	index.dxx

all:: doc

doc::
	@if [ ! -d $(DOCDIR) ]; then mkdir $(DOCDIR); fi
	@rm -f $(DOCDIR)/*.html
	cat Doxyfile.in | \
          sed -e "s/\(PROJECT_NUMBER *= \)/\1 $(VIGRA_VERSION)/" \
              -e "s@\(HTML_OUTPUT *= \).*@\1 $(DOCDIR)@" > Doxyfile
	$(DOC)
	python post.py $(DOCDIR) $(VIGRA_VERSION)
	cp $(DOCDIR)/index.html $(DOCDIR)/index.html.in
	cp ../src/examples/lenna*.gif $(DOCDIR)

latex::
	$(DOCLATEX) $(FILES)

clean::
	rm -f $(DOCDIR)/*.html $(DOCDIR)/*.gif $(DOCDIR)/*.pnm \
              $(DOCDIR)/gifs.db $(DOCDIR)/*.ps $(DOCDIR)/*.eps \
              $(DOCDIR)/dxxgifs* $(DOCDIR)/*.class $(DOCDIR)/*.in \
              $(DOCDIR)/formula.repository $(DOCDIR)/*.css \
              $(DOCDIR)/_formulas.* $(DOCDIR)/doxfont.ttf
