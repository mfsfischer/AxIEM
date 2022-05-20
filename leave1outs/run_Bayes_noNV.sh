#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/test_Bayes_classifier.py -a "$VIRUS"_testnoNV.txt -b "$VIRUS"_noNV.txt -e noNV -p BayesPrecisionRecallnoNV.txt -r BayesROC_curvesnoNV.txt -z BayesSummaryStatisticsnoNV.txt -y BayesNeighborSumValuesnoNV.txt
done
