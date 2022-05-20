#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    grep -F "$VIRUS" NeighborSumValues48.txt > "$VIRUS"_NS48.txt
    grep -F -v "$VIRUS" NeighborSumValues48.txt > "$VIRUS"_testNS48.txt
    python ../../../scripts/LinearRegression_valuesonly.py -a "$VIRUS"_testNS48.txt -b "$VIRUS"_NS48.txt -p PrecisionRecallNS48.txt -r ROC_curvesNS48.txt -z SummaryStatisticsNS48.txt -y NS48values.txt -e 48
done
