#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/LinearRegression_NeighborSum.py -a ../training_datasets/"$VIRUS"_test.csv -b ../testing_datasets/"$VIRUS"_epi.csv -e 48 -l 4 -p PrecisionRecall48.txt -r ROC_curves48.txt -z SummaryStatistics48.txt -y NeighborSumValues48.txt
done
