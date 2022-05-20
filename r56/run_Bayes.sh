#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/test_Bayes_classifier.py -a "$VIRUS"_testNS56.txt -b "$VIRUS"_NS56.txt -e 56 -p BayesPrecisionRecall56.txt -r BayesROC_curves56.txt -z BayesSummaryStatistics56.txt -y BayesNeighborSumValues56.txt
done
