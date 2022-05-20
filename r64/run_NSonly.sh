#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    grep -F "$VIRUS" NeighborSumValues64.txt > "$VIRUS"_NS64.txt
    grep -F -v "$VIRUS" NeighborSumValues64.txt > "$VIRUS"_testNS64.txt
    python ../../../scripts/LinearRegression_valuesonly.py -a "$VIRUS"_testNS64.txt -b "$VIRUS"_NS64.txt -p PrecisionRecallNS64.txt -r ROC_curvesNS64.txt -z SummaryStatisticsNS64.txt -y NS64values.txt -e 64
done
