#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/LinearRegression_NeighborSum.py -a ../training_datasets/"$VIRUS"_test.csv -b ../testing_datasets/"$VIRUS"_epi.csv -e 32 -l 4 -p PrecisionRecall32.txt -r ROC_curves32.txt -z SummaryStatistics32.txt -y NeighborSumValues32.txt
done
