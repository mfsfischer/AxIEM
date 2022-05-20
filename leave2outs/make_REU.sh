#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    grep -F "$VIRUS" NeighborSumValues_REU.txt > "$VIRUS"_REU.txt
    grep -F -v "$VIRUS" NeighborSumValues_REU.txt > "$VIRUS"_testREU.txt
done
