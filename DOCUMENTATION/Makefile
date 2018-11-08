default: network_synthesis_manual.pdf

LATEX_FILES= \
 Tex/introduction.tex \
 Tex/network_synthesis.tex \
 Tex/software_user_guide.tex \
 Tex/examples.tex \
 Tex/bibliography.tex

network_synthesis_manual.pdf: network_synthesis.tex  $(LATEX_FILES)
	pdflatex network_synthesis.tex 
	pdflatex network_synthesis.tex 
	pdflatex network_synthesis.tex 

clean:
	rm -f *.aux *.lof *.log  *.out  *.toc 

clean_all:
	rm -f *.aux *.lof *.log  *.out  *.toc *.pdf

