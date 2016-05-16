# Name of manuscript
manuscript = phdthesis
summary = exec_summ
presentation = slides

# Source files
sources = $(wildcard *.tex) \
          $(wildcard frontmatter/*.tex) \
          $(wildcard chapters/*.tex) \
          $(wildcard figures/*.tex) \
          $(wildcard algorithms/*.tex) \
          $(wildcard slide_src/*.tex) \
          references.bib

# PdfLaTeX compilation options
latexopt = -halt-on-error -file-line-error

#=================================================================
# Generate PDF of manuscript using PdfLaTeX
#=================================================================

all: $(manuscript).pdf

$(manuscript).pdf: $(sources)
	pdflatex $(latexopt) $(manuscript).tex
	bibtex $(manuscript)
	pdflatex $(latexopt) $(manuscript).tex
	pdflatex $(latexopt) $(manuscript).tex

#=================================================================
# Generate PDF of executive summary
#=================================================================

exec: $(summary).pdf

$(summary).pdf: $(sources)
	pdflatex $(latexopt) $(summary).tex
	bibtex $(summary)
	pdflatex $(latexopt) $(summary).tex
	pdflatex $(latexopt) $(summary).tex

#=================================================================
# Generate PDF of presentation
#=================================================================

slides: $(presentation).pdf

$(presentation).pdf: $(sources)
	mkdir -p presentation
	cd presentation; pdflatex --shell-escape ../slide_src/slides.tex; pdflatex --shell-escape ../slide_src/slides.tex; cp slides.pdf ..

#=================================================================
# Other
#=================================================================

clean:
	@rm -f *.aux *.bbl *.blg *.log *.out *.spl *.lof *.lot *.toc
	@rm -f chapters/*.aux

.PHONY: all clean
