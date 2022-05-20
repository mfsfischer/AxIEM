#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/randomize_predictors.py -a "$VIRUS"_testNS64.txt -b "$VIRUS"_NS64.txt -e 64 -p RandBayesPrecisionRecall64.txt -r RandBayesROC_curves64.txt -z RandBayesSummaryStatistics64.txt -y RandBayesNeighborSumValues64.txt
done
