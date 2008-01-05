ifndef docdir
    docdir=../doc/vigra
endif

         DOC = doxygen
    DOCLATEX = doc++ -u -o $(DOCDIR)/latex/vigra.tex -t
      DOCDIR = $(docdir)

FILES = \
	index.dxx

all:: doxygen post

doxygen::
	@if [ ! -d $(DOCDIR) ]; then mkdir -p $(DOCDIR); fi
	@rm -f $(DOCDIR)/*.html
	cat Doxyfile.in | \
          sed -e "s/\(PROJECT_NUMBER *= \)/\1 $(VIGRA_VERSION)/" \
              -e "s/\(vigra_version=\)/\1\"$(VIGRA_VERSION)\"/" \
              -e "s@\(HTML_OUTPUT *= \).*@\1 $(DOCDIR)@" > Doxyfile
	$(DOC)
        @rm -f doxygen*.tmp

post::
#	python post.py $(DOCDIR) $(VIGRA_VERSION)
	python makeFunctionIndex.py $(DOCDIR)
	cp $(DOCDIR)/index.html $(DOCDIR)/index.html.in
	cp ../src/examples/lenna*.gif $(DOCDIR)

latex::
	$(DOCLATEX) $(FILES)

clean::
	rm -f $(DOCDIR)/*.html $(DOCDIR)/*.gif $(DOCDIR)/*.pnm $(DOCDIR)/*.png \
              $(DOCDIR)/gifs.db $(DOCDIR)/*.ps $(DOCDIR)/*.eps \
              $(DOCDIR)/dxxgifs* $(DOCDIR)/*.class $(DOCDIR)/*.in \
              $(DOCDIR)/formula.repository $(DOCDIR)/*.css \
              $(DOCDIR)/_formulas.* $(DOCDIR)/doxfont.ttf
