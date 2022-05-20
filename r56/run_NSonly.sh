#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    grep -F "$VIRUS" NeighborSumValues56.txt > "$VIRUS"_NS56.txt
    grep -F -v "$VIRUS" NeighborSumValues56.txt > "$VIRUS"_testNS56.txt
    python ../../../scripts/LinearRegression_valuesonly.py -a "$VIRUS"_testNS56.txt -b "$VIRUS"_NS56.txt -p PrecisionRecallNS56.txt -r ROC_curvesNS56.txt -z SummaryStatisticsNS56.txt -y NS56values.txt -e 56
done
