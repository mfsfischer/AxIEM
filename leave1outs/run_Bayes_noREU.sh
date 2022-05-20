#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/test_Bayes_classifier.py -a "$VIRUS"_testnoREU.txt -b "$VIRUS"_noREU.txt -e noREU -p BayesPrecisionRecallnoREU.txt -r BayesROC_curvesnoREU.txt -z BayesSummaryStatisticsnoREU.txt -y BayesNeighborSumValuesnoREU.txt
done
