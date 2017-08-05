ps: ICC16.tex 
	latex ICC16.tex
	#biber ICC16.bcf
	#latex ICC16.tex
	dvips ICC16.dvi -o ICC16.ps
	ps2pdf ICC16.ps


pdfonly: ICC16.tex
	pdflatex ICC16.tex

pdf: ICC16.tex
	pdflatex ICC16.tex
	biber ICC16.bcf
	pdflatex ICC16.tex

clean:
	rm *.log *.aux *.bbl *.blg *.pdf *.ps *.dvi *.out *.snm *.nav *.toc *.bcf *.xml


cp_pres:
	cp -rf ICC16.tex ICC16.pdf Makefile templates/ logos ~/Dropbox/ICC2016/pres/ 	
