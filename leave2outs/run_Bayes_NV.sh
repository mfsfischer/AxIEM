#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/test_Bayes_classifier.py -a "$VIRUS"_testNV.txt -b "$VIRUS"_NV.txt -e NV -p BayesPrecisionRecallNV.txt -r BayesROC_curvesNV.txt -z BayesSummaryStatisticsNV.txt -y BayesNeighborSumValuesNV.txt
done
