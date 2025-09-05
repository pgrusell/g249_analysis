TARGET := report
ENGINE ?= pdf        
BIB    ?= bibtex     

TEXSRC := $(TARGET).tex

LATEXMK_BIBFLAG := $(if $(filter $(BIB),biber),-use-biber,-bibtex)

.PHONY: all clean distclean

all: $(TARGET).pdf

$(TARGET).pdf: $(TEXSRC)
	latexmk -$(ENGINE) -interaction=nonstopmode -halt-on-error \
		-file-line-error -synctex=1 $(LATEXMK_BIBFLAG) $(TEXSRC)

clean:
	latexmk -C
	rm -f *.bbl *.blg

distclean: clean
	rm -f $(TARGET).pdf