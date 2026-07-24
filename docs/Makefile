# Minimal makefile for the bilingual sqc documentation site.
# `make html` builds English at build/html and Chinese at build/html/zh.
SPHINXOPTS  ?= -W --keep-going
SPHINXBUILD ?= sphinx-build
SOURCEDIR    = source
BUILDDIR     = build

.PHONY: help html en zh clean

help:
	@echo "make html   - build both languages (en at root, zh at /zh/)"
	@echo "make en      - build English only"
	@echo "make zh      - build Chinese only"
	@echo "make clean   - remove build output"

# The tutorial notebook lives at the repo root (single source of truth). Copy it
# into each language tree before building so nbsphinx can pick it up. The copies
# are gitignored.
NOTEBOOK = ../Simulation_sqc.ipynb

notebooks:
	cp "$(NOTEBOOK)" "$(SOURCEDIR)/en/Simulation_sqc.ipynb"
	cp "$(NOTEBOOK)" "$(SOURCEDIR)/zh/Simulation_sqc.ipynb"

en: notebooks
	$(SPHINXBUILD) $(SPHINXOPTS) -b html "$(SOURCEDIR)/en" "$(BUILDDIR)/html"

zh: notebooks
	$(SPHINXBUILD) $(SPHINXOPTS) -b html "$(SOURCEDIR)/zh" "$(BUILDDIR)/html/zh"

html: en zh

clean:
	rm -rf "$(BUILDDIR)"
	rm -f "$(SOURCEDIR)/en/Simulation_sqc.ipynb" "$(SOURCEDIR)/zh/Simulation_sqc.ipynb"
