#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/test_Bayes_classifier.py -a "$VIRUS"_testNS16.txt -b "$VIRUS"_NS16.txt -e 16 -p BayesPrecisionRecall16.txt -r BayesROC_curves16.txt -z BayesSummaryStatistics16.txt -y BayesNeighborSumValues16.txt
done
