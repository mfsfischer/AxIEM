#!/bin/bash

for VIRUS in $(cat virus.list)

do
	python LinearRegression_NeighborSum.py -a "$VIRUS"_test_epi.csv -b "$VIRUS"_epi.csv -e 16 -p PrecisionRecall.txt -r ROC_curves.txt -z SummaryStatistics.txt -y NeighborSumValues.txt
done
