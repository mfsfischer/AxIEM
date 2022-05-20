#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/test_Bayes_classifier.py -a "$VIRUS"_testCP.txt -b "$VIRUS"_CP.txt -e CP -p BayesPrecisionRecallCP.txt -r BayesROC_curvesCP.txt -z BayesSummaryStatisticsCP.txt -y BayesNeighborSumValuesCP.txt
done
