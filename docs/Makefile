# You can set these variables from the command line.
SPHINXOPTS    = -n -j 8
PAPER         =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = Finesse
SOURCEDIR     = source
BUILDDIR      = build
EXECDIR       = $(BUILDDIR)/jupyter_execute
RESLATE       = reslate
PKGDIR        = ../src/finesse/
SRCDIRS       = $(shell find ../src/finesse/ -type d)
SRCFILES      = $(shell find ../src/finesse/ -type f -name '*')

# internal
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS)

# prevents warnings about frozen modules clogging up the output
export PYDEVD_DISABLE_FILE_VALIDATION = 1

# make sure this is first so "make" = "make help"
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help clean apipages html latex latexpdf man epub text htmlfast htmlsingleproc

clean:
	-rm -f .apipages
	-rm -rf $(BUILDDIR)/* ./html
	find $(SOURCEDIR)/api -type f -name "*.rst" ! -name "index.rst" -delete
	find $(SOURCEDIR)/api -type d -empty -delete

IPYNB_SRC = $(shell find $(SOURCEDIR) -name "*.ipynb" -not -path '*/.*')
IPYNB_TRG = $(subst $(SOURCEDIR),$(EXECDIR),$(IPYNB_SRC:.ipynb=.py))

ipynb2python: $(IPYNB_TRG)

$(EXECDIR)/%.py : $(SOURCEDIR)/%.ipynb
	jupyter-nbconvert --to=python --no-prompt --output-dir=$(dir $@) $<

.apipages: $(SRCFILES)
	$(RESLATE) --clean
	$(RESLATE) $(PKGDIR) -v

apipages: .apipages

html_build: apipages ipynb2python
	$(SPHINXBUILD) -M html $(SOURCEDIR) $(BUILDDIR) $(ALLSPHINXOPTS) -t html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

html: html_build
	touch .apipages

htmlfast_build: ipynb2python
	$(SPHINXBUILD) -M html $(SOURCEDIR) $(BUILDDIR) $(ALLSPHINXOPTS) -t htmlfast
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

htmlfast: htmlfast_build
	touch .apipages

htmlsingleproc_build: apipages ipynb2python
	$(SPHINXBUILD) -M html $(SOURCEDIR) $(BUILDDIR) $(sed "s/-j auto//g" <<< $(ALLSPHINXOPTS)) -t htmlsingleproc
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

htmlsingleproc: htmlsingleproc_build

latex_build:
	$(SPHINXBUILD) -M latex $(SOURCEDIR) $(BUILDDIR) $(ALLSPHINXOPTS) -t latex
	@echo
	@echo "Build finished. The LaTeX files are in $(BUILDDIR)/latex."

latex: latex_build

latexpdf_build:
	$(SPHINXBUILD) -M latex $(SOURCEDIR) $(BUILDDIR) $(ALLSPHINXOPTS) -t latex
	@echo
	@echo "Build finished. The LaTeX files are in $(BUILDDIR)/latex."
	@echo "Running make on the LaTeX files in $(BUILDDIR)/latex..."
	$(MAKE) -C $(BUILDDIR)/latex all-pdf

latexpdf: latexpdf_build

man_build: apipages
	$(SPHINXBUILD) -M man $(SOURCEDIR) $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "Build finished. The manual pages are in $(BUILDDIR)/man."

man: man_build
	touch .apipages

epub_build: apipages
	$(SPHINXBUILD) -M epub $(ALLSPHINXOPTS) $(BUILDDIR)/epub
	@echo
	@echo "Build finished. The epub file is in $(BUILDDIR)/epub."

epub: epub_build
	touch .apipages

text_build: apipages
	$(SPHINXBUILD) -M text $(ALLSPHINXOPTS) $(BUILDDIR)/text
	@echo
	@echo "Build finished. The text files are in $(BUILDDIR)/text."

text: text_build
	touch .apipages
