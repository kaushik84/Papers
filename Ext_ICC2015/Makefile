ps_1C: paper_1C.tex refs.bib
	latex paper_1C.tex
	bibtex paper_1C
	latex paper_1C.tex
	dvips paper_1C.dvi -o paper_1C.ps
	ps2pdf paper_1C.ps

pdf_1C: paper_1C.tex refs.bib
	pdflatex paper_1C.tex
	bibtex paper_1C
	pdflatex paper_1C.tex

# Only for embedding pdf fonts
pdf_EF_1C: paper_1C.tex refs.bib 
	pdflatex paper_1C.tex
	bibtex paper_1C
	pdflatex paper_1C.tex
	# convert to Postscript
	pdftops paper_1C.pdf paper_1C.ps
	# reconvert to PDF but enforce font embedding (note, hyperlinks will not survive)
	ps2pdf14 -dPDFSETTINGS=/prepress -dEmbedAllFonts=true paper_2C.ps paper_2C.pdf
	# verify
	pdffonts paper_1C.pdf

clean_1C:
	rm *.log *.aux *.bbl *.blg paper_1C.pdf *.ps *.dvi *.out

cp_complete_1C:
	cp -rf ../Ext_CrownCom2015/ ~/Dropbox 	


cp_paper_1C:
	cp -rf paper_1C.tex paper_1C.pdf refs.bib figures/ ~/Dropbox/Ext_ICC2015/ 

	
ps_2C: paper_2C.tex refs.bib
	latex paper_2C.tex
	bibtex paper_2C
	latex paper_2C.tex
	dvips paper_2C.dvi -o paper_2C.ps
	ps2pdf paper_2C.ps

pdf_2C: paper_2C.tex refs.bib
	pdflatex paper_2C.tex
	bibtex paper_2C
	pdflatex paper_2C.tex

# Only for embedding pdf fonts
pdf_EF_2C: paper_2C.tex refs.bib 
	pdflatex paper_2C.tex
	bibtex paper_2C
	pdflatex paper_2C.tex
	# convert to Postscript
	pdftops paper_2C.pdf paper_2C.ps
	# reconvert to PDF but enforce font embedding (note, hyperlinks will not survive)
	ps2pdf14 -dPDFSETTINGS=/prepress -dEmbedAllFonts=true paper_2C.ps paper_2C.pdf
	# verify
	pdffonts paper_2C.pdf

clean_2C:
	rm *.log *.aux *.bbl *.blg paper_2C.pdf *.ps *.dvi *.out

cp_complete_2C:
	cp -rf ../Ext_CrownCom2015/ ~/Dropbox 	


cp_paper_2C:
	cp -rf paper_2C.tex paper_2C.pdf refs.bib figures/ ~/Dropbox/Ext_ICC2015/ 	
