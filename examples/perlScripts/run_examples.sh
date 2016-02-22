#!/bin/bash

echo "===================================================================="

echo "ltcstool.pl:"
echo "  start ltcstool.pl and save results in ltcs_*.csv, LPs in ltcs_*.lp"
read -n 1 -p "(press any key)"

perl ../../scripts/ltcstool.pl -s model.xml -r rfile.txt -e modes.example -c conc.csv -g dfg.txt -t 4 -l ltcs_ -f ltcs_ -d 

echo "--------------------------------------------------------------------"

echo "ltcstool_hotstart.pl:"
echo "  start ltcstool_hotstart.pl with temp_4.lp and intervall 10 and"
echo "  save results in hs_*.csv, LPs in hs_*.lp"
read -n 1 -p "(press any key)"

perl ../../scripts/ltcstool_hotstart.pl -i ltcs_4.lp -n 10 -l hs_ -f hs_ -t 4 -d 

echo "--------------------------------------------------------------------"

echo "ltcstool_clever.pl:"
echo "  start ltcstool_clever.pl using timeout 0.5 and save results in"
echo "  clever_*.csv, LPs in clever_*.lp"
read -n 1 -p "(press any key)"

perl ../../scripts/ltcstool_clever.pl -s model.xml -r rfile.txt -e modes.example -c conc.csv -g dfg.txt -x 0.5 -l clever_ -f clever_ -d 

echo "--------------------------------------------------------------------"

echo "ltcstool_hotstart_clever.pl:"
echo "  start ltcstool_hotstart_clever.pl with temp_4.lp and timeout 0.5"
echo "  and save results in hs_clever_*.csv, LPs in hs_clever_*.lp"
read -n 1 -p "(press any key)"

perl ../../scripts/ltcstool_hotstart_clever.pl -i ltcs_4.lp -x 0.5 -l hs_clever_ -f hs_clever_ -d 

echo "--------------------------------------------------------------------"

echo "ltcstool.pl with yields:"
echo "  start ltcstool.pl using yield information and save results in"
echo "  ltcs_y_*.csv, LPs in ltcs_y_*.lp"
read -n 1 -p "(press any key)"

perl ../../scripts/ltcstool.pl -s model.xml -r rfile.txt -e modes.example -c conc.csv -g dfg.txt -t 4 -y yields.txt -l ltcs_y_ -f ltcs_y_ -d 

echo "--------------------------------------------------------------------"

echo "ltcstool_clever.pl with yields:"
echo "  start ltcstool_clever.pl using using yield information and save"
echo "  results in clever_y_*.csv, LPs in clever_y_*.lp"
read -n 1 -p "(press any key)"

perl ../../scripts/ltcstool_clever.pl -s model.xml -r rfile.txt -e modes.example -c conc.csv -g dfg.txt -y yields.txt -l clever_y_ -f clever_y_ -d 

echo "===================================================================="

