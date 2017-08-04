all:	draft.tex refs.bib
	latex draft.tex
	bibtex draft
	latex draft.tex
	dvips draft.dvi -o draft.ps
	ps2pdf draft.ps

pdf:    draft.tex refs.bib
	pdflatex draft.tex
	bibtex draft
	pdflatex draft.tex

# Only for embedding pdf fonts
pdf_EF: draft.tex refs.bib 
	pdflatex draft.tex
	bibtex draft
	pdflatex draft.tex
	# convert to Postscript
	pdftops draft.pdf draft.ps
	# reconvert to PDF but enforce font embedding (note, hyperlinks will not survive)
	ps2pdf14 -dPDFSETTINGS=/prepress -dEmbedAllFonts=true draft.ps draft.pdf
	# verify
	pdffonts draft.pdf

clean:
	rm *.log *.aux *.bbl *.blg draft.pdf *.ps *.dvi *.out



