#!/bin/bash
dlats="(/15,20,25/)"
dlatn="(/40,45,50/)"
dlonl="(/20,105,250/)"
dlonr="(/60,175,320/)"

#ncl -nQ lats=$dlats latn=$dlatn lonl=$dlonl lonr=$dlonr 2002-calc_eof_6kinds.ncl
#ncl 2002-draw_eof_regression_horizontal.ncl
#ncl 2002-draw_eof_regression_vertical.ncl
ncl -nQ lats=$dlats latn=$dlatn lonl=$dlonl lonr=$dlonr 2102-calc_corr_project_dzdt.ncl
python 2112-draw_heatmap_corr_project_dzdt.py

