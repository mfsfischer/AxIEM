#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/randomize_predictors.py -a "$VIRUS"_testNS16.txt -b "$VIRUS"_NS16.txt -e 16 -p RandBayesPrecisionRecall16.txt -r RandBayesROC_curves16.txt -z RandBayesSummaryStatistics16.txt -y RandBayesNeighborSumValues16.txt
done
