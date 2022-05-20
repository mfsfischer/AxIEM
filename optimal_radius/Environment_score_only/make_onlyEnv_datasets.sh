#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    cat "$VIRUS"_NSLM.txt | awk '{print $1,$2,$6,$6,$7,$8,$9}' > "$VIRUS"_onlyLM.txt
    cat "$VIRUS"_testNSLM.txt | awk '{print $1,$2,$6,$6,$7,$8,$9}' > "$VIRUS"_testonlyLM.txt
    python ../../../scripts/test_Bayes_classifier.py -a "$VIRUS"_testonlyLM.txt -b "$VIRUS"_onlyLM.txt -e LM -p BayesPrecisionRecallLM_Envir.txt -r BayesROC_curvesLM_Envir.txt -z BayesSummaryStatisticsLM_Envir.txt -y BayesNeighborSumValuesLM_Envir.txt
done
