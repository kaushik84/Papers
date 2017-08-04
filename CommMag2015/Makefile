ps: paper.tex refs.bib
	latex paper.tex
	bibtex paper
	latex paper.tex
	dvips paper.dvi -o paper.ps
	ps2pdf paper.ps
	

ps_EF: paper.tex refs.bib
	latex paper.tex
	bibtex paper
	latex paper.tex
	dvips paper.dvi -o paper.ps
	ps2pdf14 -dPDFSETTINGS=/prepress -dEmbedAllFonts=true paper.ps paper.pdf

pdf: paper.tex refs.bib
	pdflatex paper.tex
	bibtex paper
	pdflatex paper.tex

# Only for embedding pdf fonts
pdf_EF: paper.tex refs.bib 
	pdflatex paper.tex
	bibtex paper
	pdflatex paper.tex
	# convert to Postscript
	pdftops paper.pdf paper.ps
	# reconvert to PDF but enforce font embedding (note, hyperlinks will not survive)
	ps2pdf14 -dPDFSETTINGS=/prepress -dEmbedAllFonts=true paper.ps paper.pdf
	# verify
	pdffonts paper.pdf

clean:
	rm *.log *.aux *.bbl *.blg paper.pdf *.ps *.dvi *.out

cp_complete:
	cp -rf ../CommMag2015/ ~/Dropbox 	


cp_paper:
	cp -rf paper.tex paper.pdf refs.bib figures/  ~/Dropbox/CommMag/ 	
