all: chapter.tex refs.bib
	latex chapter.tex
	bibtex chapter 
	latex chapter.tex
	latex chapter.tex
	dvips chapter.dvi -o chapter.ps
	ps2pdf chapter.ps
clean: 
	rm *.log *.aux *.bbl *.blg chapter.pdf *.ps *.dvi *.out *.tcp *.idx *.toc *.tps

