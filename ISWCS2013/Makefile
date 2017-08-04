ps: paper.tex refs.bib
	latex paper.tex
	bibtex paper
	latex paper.tex
	dvips paper.dvi -o paper.ps
	ps2pdf paper.ps

pdf: paper.tex refs.bib
	pdflatex paper.tex
	bibtex paper
	pdflatex paper.tex

clean:
	rm *.log *.aux *.bbl *.blg *.pdf *.ps *.dvi *.out

