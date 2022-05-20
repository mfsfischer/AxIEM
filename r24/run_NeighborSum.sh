#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/LinearRegression_NeighborSum.py -a ../training_datasets/"$VIRUS"_test.csv -b ../testing_datasets/"$VIRUS"_epi.csv -e 24 -l 4 -p PrecisionRecall24.txt -r ROC_curves24.txt -z SummaryStatistics24.txt -y NeighborSumValues24.txt
done
