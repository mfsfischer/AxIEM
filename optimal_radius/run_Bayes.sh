#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    python ../../../scripts/test_Bayes_classifier.py -a "$VIRUS"_testNSLM.txt -b "$VIRUS"_NSLM.txt -e LM -p BayesPrecisionRecallLM.txt -r BayesROC_curvesLM.txt -z BayesSummaryStatisticsLM.txt -y BayesNeighborSumValuesLM.txt
done
