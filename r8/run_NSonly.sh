#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    grep -F "$VIRUS" NeighborSumValues8.txt > "$VIRUS"_NS8.txt
    grep -F -v "$VIRUS" NeighborSumValues8.txt > "$VIRUS"_testNS8.txt
    python ../../../scripts/LinearRegression_valuesonly.py -a "$VIRUS"_testNS8.txt -b "$VIRUS"_NS8.txt -p PrecisionRecallNS8.txt -r ROC_curvesNS8.txt -z SummaryStatisticsNS8.txt -y NS8values.txt -e 8
done
