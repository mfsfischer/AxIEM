#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/LinearRegression_NeighborSum.py -a ../training_datasets/"$VIRUS"_test.csv -b ../testing_datasets/"$VIRUS"_epi.csv -e 16 -l 4 -p PrecisionRecall16.txt -r ROC_curves16.txt -z SummaryStatistics16.txt -y NeighborSumValues16.txt
done
