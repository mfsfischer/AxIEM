#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/test_Bayes_classifier.py -a "$VIRUS"_testNS48.txt -b "$VIRUS"_NS48.txt -e 48 -p BayesPrecisionRecall48.txt -r BayesROC_curves48.txt -z BayesSummaryStatistics48.txt -y BayesNeighborSumValues48.txt
done
