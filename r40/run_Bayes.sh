#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/test_Bayes_classifier.py -a "$VIRUS"_testNS40.txt -b "$VIRUS"_NS40.txt -e 40 -p BayesPrecisionRecall40.txt -r BayesROC_curves40.txt -z BayesSummaryStatistics40.txt -y BayesNeighborSumValues40.txt
done
