#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/randomize_predictors.py -a "$VIRUS"_testNS40.txt -b "$VIRUS"_NS40.txt -e 40 -p RandBayesPrecisionRecall40.txt -r RandBayesROC_curves40.txt -z RandBayesSummaryStatistics40.txt -y RandBayesNeighborSumValues40.txt
done
