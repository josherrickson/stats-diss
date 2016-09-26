MAIN=errickson_thesis

default:
	R -e "knitr::knit('$(MAIN).Rnw')"
	latexmk -pdf -pv -pdflatex="pdflatex %O %S" $(MAIN).tex

bg:
	R -e "knitr::knit('$(MAIN).Rnw')"
	latexmk -pdf -pvc -pdflatex="pdflatex %O %S" $(MAIN).tex

fresh: clean default

open:
	open $(MAIN).pdf

dependencies:
	@R -e "file <- '$(MAIN).Rnw'; source('dependencies.R')"

clean:
	@rm -f $(MAIN).{tex,bbl,aux,log,lof,lot,toc,blg,out,fdb_latexmk,fls,loa,loap,synctex.gz}
	@rm -f $(MAIN)-concordance.tex .Rhistory
	@rm -rf figure/

clean-all: clean
	@rm -f $(MAIN).pdf
