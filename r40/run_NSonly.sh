#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    grep -F "$VIRUS" NeighborSumValues40.txt > "$VIRUS"_NS40.txt
    grep -F -v "$VIRUS" NeighborSumValues40.txt > "$VIRUS"_testNS40.txt
    python ../../../scripts/LinearRegression_valuesonly.py -a "$VIRUS"_testNS40.txt -b "$VIRUS"_NS40.txt -p PrecisionRecallNS40.txt -r ROC_curvesNS40.txt -z SummaryStatisticsNS40.txt -y NS40values.txt -e 40
done
