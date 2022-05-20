#!/bin/bash

VIRUS=(EBOV H3 H7 RSV SARS1 SARS2 HIV)
RADIUS=(25.6 20.1 36.1 34.9 58.9 55.0 34.9)

for ((i=0; i<${#VIRUS[@]}; i++)); do
    python LinearRegression_NeighborSum.py -a ${VIRUS[i]}_test_noNV.txt -b ${VIRUS[i]}_noNV.txt -e ${RADIUS[i]} -l 4 -p PrecisionRecall_noNV.txt -r ROC_curves_noNV.txt -z SummaryStatistics_noNV.txt -y NeighborSumValues_noNV.txt
done
