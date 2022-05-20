#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/LinearRegression_NeighborSum.py -a ../training_datasets/"$VIRUS"_test.csv -b ../testing_datasets/"$VIRUS"_epi.csv -e 8 -l 4 -p PrecisionRecall8.txt -r ROC_curves8.txt -z SummaryStatistics8.txt -y NeighborSumValues8.txt
done
