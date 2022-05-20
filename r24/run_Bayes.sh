#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/test_Bayes_classifier.py -a "$VIRUS"_testNS24.txt -b "$VIRUS"_NS24.txt -e 24 -p BayesPrecisionRecall24.txt -r BayesROC_curves24.txt -z BayesSummaryStatistics24.txt -y BayesNeighborSumValues24.txt
done
