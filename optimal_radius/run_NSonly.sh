#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    grep -F "$VIRUS" NeighborSumValuesLM.txt > "$VIRUS"_NSLM.txt
    grep -F -v "$VIRUS" NeighborSumValuesLM.txt > "$VIRUS"_testNSLM.txt
    python ../../../scripts/LinearRegression_valuesonly.py -a "$VIRUS"_testNSLM.txt -b "$VIRUS"_NSLM.txt -p PrecisionRecallNSLM.txt -r ROC_curvesNSLM.txt -z SummaryStatisticsNSLM.txt -y NSLMvalues.txt -e LM
done
