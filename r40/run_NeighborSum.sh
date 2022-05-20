#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/LinearRegression_NeighborSum.py -a ../training_datasets/"$VIRUS"_test.csv -b ../testing_datasets/"$VIRUS"_epi.csv -e 40 -l 4 -p PrecisionRecall40.txt -r ROC_curves40.txt -z SummaryStatistics40.txt -y NeighborSumValues40.txt
done
