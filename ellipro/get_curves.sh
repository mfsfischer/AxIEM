#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/get_PR-ROC_curve_stats.py -d "$VIRUS"_ellipro.txt -e Ellipro -p PrecisionRecall_ellipro.txt -r ROC_curves_ellipro.txt 
done
