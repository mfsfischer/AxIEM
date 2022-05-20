#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    grep -F "$VIRUS" NeighborSumValues_NV.txt > "$VIRUS"_NV.txt
    grep -F -v "$VIRUS" NeighborSumValues_NV.txt > "$VIRUS"_testNV.txt
done
