#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/test_Bayes_classifier.py -a "$VIRUS"_testNS8.txt -b "$VIRUS"_NS8.txt -e 8 -p BayesPrecisionRecall8.txt -r BayesROC_curves8.txt -z BayesSummaryStatistics8.txt -y BayesNeighborSumValues8.txt
done
