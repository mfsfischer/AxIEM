#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    grep -F "$VIRUS" NeighborSumValues32.txt > "$VIRUS"_NS32.txt
    grep -F -v "$VIRUS" NeighborSumValues32.txt > "$VIRUS"_testNS32.txt
    python ../../../scripts/LinearRegression_valuesonly.py -a "$VIRUS"_testNS32.txt -b "$VIRUS"_NS32.txt -p PrecisionRecallNS32.txt -r ROC_curvesNS32.txt -z SummaryStatisticsNS32.txt -y NS32values.txt -e 32
done
