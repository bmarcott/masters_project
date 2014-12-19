all: project-report.pdf

project-report.pdf: project-report.tex
	pdflatex project-report.tex
	bibtex project-report
	pdflatex project-report.tex
	pdflatex project-report.tex

.PHONY: clean cleanest

clean:
	rm -f *.aux *.log *.bbl *.blg

cleanest: clean
	rm -f project-report.pdf


