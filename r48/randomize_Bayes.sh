#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/randomize_predictors.py -a "$VIRUS"_testNS48.txt -b "$VIRUS"_NS48.txt -e 48 -p RandBayesPrecisionRecall48.txt -r RandBayesROC_curves48.txt -z RandBayesSummaryStatistics48.txt -y RandBayesNeighborSumValues48.txt
done
