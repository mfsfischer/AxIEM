#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/randomize_predictors.py -a "$VIRUS"_testNS8.txt -b "$VIRUS"_NS8.txt -e 8 -p RandBayesPrecisionRecall8.txt -r RandBayesROC_curves8.txt -z RandBayesSummaryStatistics8.txt -y RandBayesNeighborSumValues8.txt
done
