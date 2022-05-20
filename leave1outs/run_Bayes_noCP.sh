#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/test_Bayes_classifier.py -a "$VIRUS"_testnoCP.txt -b "$VIRUS"_noCP.txt -e noCP -p BayesPrecisionRecallnoCP.txt -r BayesROC_curvesnoCP.txt -z BayesSummaryStatisticsnoCP.txt -y BayesNeighborSumValuesnoCP.txt
done
