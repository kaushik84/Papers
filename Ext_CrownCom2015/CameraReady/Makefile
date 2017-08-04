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
	cp -rf paper_2C.tex paper_2C.pdf refs.bib figures/ SourceFiles ~/Dropbox/Ext_CrownCom2015/ 	



cp_resp:
	cp -rf ./response\ letter/* ~/Dropbox/Ext_CrownCom2015/response\ letter/	
