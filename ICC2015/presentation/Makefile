ps: ICC15.tex 
	latex ICC15.tex
	#biber ICC15.bcf
	#latex ICC15.tex
	dvips ICC15.dvi -o ICC15.ps
	ps2pdf ICC15.ps


pdfonly: ICC15.tex
	pdflatex ICC15.tex

pdf: ICC15.tex
	pdflatex ICC15.tex
	biber ICC15.bcf
	pdflatex ICC15.tex

clean:
	rm *.log *.aux *.bbl *.blg *.pdf *.ps *.dvi *.out *.snm *.nav *.toc *.bcf *.xml


cp_pres:
	cp -rf ICC15.tex ICC15.pdf Makefile templates/ logos ~/Dropbox/ICC2015/presentation/ 	
