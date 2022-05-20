#!/bin/bash

VIRUS=(EBOV H3 H7 RSV SARS1 SARS2 HIV)
RADIUS=(25.6 20.1 36.1 34.9 58.9 55.0 34.9)

for ((i=0; i<${#VIRUS[@]}; i++)); do
    python ../../../scripts/LinearRegression_NeighborSum.py -a ../training_datasets/${VIRUS[i]}_test.csv -b ../testing_datasets/${VIRUS[i]}_epi.csv -e ${RADIUS[i]} -l 4 -p PrecisionRecallLM.txt -r ROC_curvesLM.txt -z SummaryStatisticsLM.txt -y NeighborSumValuesLM.txt
done
