#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    grep -F "$VIRUS" NeighborSumValues_noREU.txt > "$VIRUS"_noREU.txt
    grep -F -v "$VIRUS" NeighborSumValues_noREU.txt > "$VIRUS"_testnoREU.txt
done
