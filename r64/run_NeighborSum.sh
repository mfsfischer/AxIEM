#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/LinearRegression_NeighborSum.py -a ../training_datasets/"$VIRUS"_test.csv -b ../testing_datasets/"$VIRUS"_epi.csv -e 64 -l 4 -p PrecisionRecall64.txt -r ROC_curves64.txt -z SummaryStatistics64.txt -y NeighborSumValues64.txt
done
