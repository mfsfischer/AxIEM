#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    grep -F "$VIRUS" NeighborSumValues16.txt > "$VIRUS"_NS16.txt
    grep -F -v "$VIRUS" NeighborSumValues16.txt > "$VIRUS"_testNS16.txt
    python ../../../scripts/LinearRegression_valuesonly.py -a "$VIRUS"_testNS16.txt -b "$VIRUS"_NS16.txt -p PrecisionRecallNS16.txt -r ROC_curvesNS16.txt -z SummaryStatisticsNS16.txt -y NS16values.txt -e 16
done
