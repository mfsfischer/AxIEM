#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/randomize_predictors.py -a "$VIRUS"_testNS24.txt -b "$VIRUS"_NS24.txt -e 24 -p RandBayesPrecisionRecall24.txt -r RandBayesROC_curves24.txt -z RandBayesSummaryStatistics24.txt -y RandBayesNeighborSumValues24.txt
done
