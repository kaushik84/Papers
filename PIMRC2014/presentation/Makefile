ps: PIMRC14.tex refs.bib
	latex PIMRC14.tex
	bibtex PIMRC14
	latex PIMRC14.tex
	dvips PIMRC14.dvi -o paper.ps
	ps2pdf PIMRC14.ps

pdfonly: PIMRC14.tex
	pdflatex PIMRC14.tex

pdf: PIMRC14.tex
	pdflatex PIMRC14.tex
	biber PIMRC14.bcf
	pdflatex PIMRC14.tex

clean:
	rm *.log *.aux *.bbl *.blg *.pdf *.ps *.dvi *.out *.snm *.nav *.toc *.bcf *.xml


