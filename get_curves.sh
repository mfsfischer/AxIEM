#!/bin/bash

for VIRUS in $(cat virus.list)

do
    python ../scripts/get_PR-ROC_curve_stats.py -d "$VIRUS"_NS32.txt -e 32 -p PR32.txt -r ROC32.txt 
    python ../scripts/get_PR-ROC_curve_stats.py -d "$VIRUS"_NS32.txt -e 32 -p PR32-NSonly.txt -r ROC32-NSonly.txt 
done
