ps: VTC16.tex 
	latex VTC16.tex
	#biber VTC16.bcf
	#latex VTC16.tex
	dvips VTC16.dvi -o VTC16.ps
	ps2pdf VTC16.ps


pdfonly: VTC16.tex
	pdflatex VTC16.tex

pdf: VTC16.tex
	pdflatex VTC16.tex
	biber VTC16.bcf
	pdflatex VTC16.tex

clean:
	rm *.log *.aux *.bbl *.blg *.pdf *.ps *.dvi *.out *.snm *.nav *.toc *.bcf *.xml


cp_pres:
	cp -rf VTC16.tex VTC16.pdf Makefile templates/ logos ~/Dropbox/VTC2016/pres/ 	
