#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    grep -F "$VIRUS" NeighborSumValues_noNV.txt > "$VIRUS"_noNV.txt
    grep -F -v "$VIRUS" NeighborSumValues_noNV.txt > "$VIRUS"_testnoNV.txt
done
