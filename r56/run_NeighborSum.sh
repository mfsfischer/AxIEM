#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/LinearRegression_NeighborSum.py -a ../training_datasets/"$VIRUS"_test.csv -b ../testing_datasets/"$VIRUS"_epi.csv -e 56 -l 4 -p PrecisionRecall56.txt -r ROC_curves56.txt -z SummaryStatistics56.txt -y NeighborSumValues56.txt
done
