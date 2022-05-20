#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/test_Bayes_classifier.py -a "$VIRUS"_testNS32.txt -b "$VIRUS"_NS32.txt -e 32 -p BayesPrecisionRecall32.txt -r BayesROC_curves32.txt -z BayesSummaryStatistics32.txt -y BayesNeighborSumValues32.txt
done
