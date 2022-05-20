#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/test_Bayes_classifier.py -a "$VIRUS"_testREU.txt -b "$VIRUS"_REU.txt -e REU -p BayesPrecisionRecallREU.txt -r BayesROC_curvesREU.txt -z BayesSummaryStatisticsREU.txt -y BayesNeighborSumValuesREU.txt
done
