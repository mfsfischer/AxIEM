#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    grep -F "$VIRUS" NeighborSumValues_noCP.txt > "$VIRUS"_noCP.txt
    grep -F -v "$VIRUS" NeighborSumValues_noCP.txt > "$VIRUS"_testnoCP.txt
done
