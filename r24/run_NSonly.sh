#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    grep -F "$VIRUS" NeighborSumValues24.txt > "$VIRUS"_NS24.txt
    grep -F -v "$VIRUS" NeighborSumValues24.txt > "$VIRUS"_testNS24.txt
    python ../../../scripts/LinearRegression_valuesonly.py -a "$VIRUS"_testNS24.txt -b "$VIRUS"_NS24.txt -p PrecisionRecallNS24.txt -r ROC_curvesNS24.txt -z SummaryStatisticsNS24.txt -y NS24values.txt -e 24
done
