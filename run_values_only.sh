#!/bin/bash

for VIRUS in $(cat virus.list)

do
	python LinearRegression_valuesonly.py -a "$VIRUS"_test_epi.csv -b "$VIRUS"_epi.csv -p PrecisionRecall.txt -r ROC_curves.txt -z SummaryStatistics.txt -y Predicted_values.txt
done
