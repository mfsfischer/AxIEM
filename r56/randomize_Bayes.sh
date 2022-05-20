#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/randomize_predictors.py -a "$VIRUS"_testNS56.txt -b "$VIRUS"_NS56.txt -e 56 -p RandBayesPrecisionRecall56.txt -r RandBayesROC_curves56.txt -z RandBayesSummaryStatistics56.txt -y RandBayesNeighborSumValues56.txt
done
