#!/bin/bash
# Proper header for a Bash script
# args: subject dispersion isGeorgeson #bootIters pmfType [sfFlag]
# pmfType - 1 for weibull, 2 for normCDF

python pmf_plot.py 2 1 1 100 2 sf15
python pmf_plot.py 2 1 1 100 2 sf30
python pmf_plot.py 2 1 1 100 2 sf60

python pmf_plot.py 2 1 1 100 1 sf15
python pmf_plot.py 2 1 1 100 1 sf30
python pmf_plot.py 2 1 1 100 1 sf60

python pmf_plot.py 1 1 1 100 2 sf15
python pmf_plot.py 1 1 1 100 2 sf30
python pmf_plot.py 1 1 1 100 2 sf60

python pmf_plot.py 1 1 1 100 1 sf15
python pmf_plot.py 1 1 1 100 1 sf30
python pmf_plot.py 1 1 1 100 1 sf60


#python pmf_plot.py 1 1 1 100 sf15
#python pmf_plot.py 1 1 1 100 sf30
#python pmf_plot.py 1 1 1 100 sf60
#python pmf_plot.py 1 3 1 100 sf15
#python pmf_plot.py 1 3 1 100 sf30
#python pmf_plot.py 2 1 0 100 
#python pmf_plot.py 2 1 1 100 
#python pmf_plot.py 2 1 1 100 sf15
#python pmf_plot.py 2 2 1 100 
#python pmf_plot.py 2 3 1 100 