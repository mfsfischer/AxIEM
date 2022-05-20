#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/test_Bayes_classifier.py -a "$VIRUS"_testNS64.txt -b "$VIRUS"_NS64.txt -e 64 -p BayesPrecisionRecall64.txt -r BayesROC_curves64.txt -z BayesSummaryStatistics64.txt -y BayesNeighborSumValues64.txt
done
