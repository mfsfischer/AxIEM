#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/get_PR-ROC_curve_stats.py -d "$VIRUS"_discotope.txt -e Discotope -p PRdiscotope.txt -r ROCdiscotope.txt 
done
