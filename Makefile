.PHONY: analysis figures clean

analysis:
	Rscript scripts/run_analysis.R

figures:
	python3 scripts/generate_report_figures.py

clean:
	rm -f results/figures/*_simulated.png
	rm -f results/tables/*_simulated.csv
