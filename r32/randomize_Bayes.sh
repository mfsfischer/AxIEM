#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/randomize_predictors.py -a "$VIRUS"_testNS32.txt -b "$VIRUS"_NS32.txt -e 32 -p RandBayesPrecisionRecall32.txt -r RandBayesROC_curves32.txt -z RandBayesSummaryStatistics32.txt -y RandBayesNeighborSumValues32.txt
done
