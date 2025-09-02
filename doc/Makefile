TARGET := report
TEXSRC := $(TARGET).tex

all: pdf

pdf: $(TARGET).pdf

$(TARGET).pdf: $(TEXSRC)
	pdflatex -interaction=nonstopmode -halt-on-error $(TEXSRC)

clean:
	rm -f *.aux *.log *.bbl *.blg *.toc *.out *.synctex.gz *.fls *.fdb_latexmk

distclean: clean
	rm -f $(TARGET).pdf