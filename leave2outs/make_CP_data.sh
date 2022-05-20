#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    grep -F "$VIRUS" NeighborSumValues_CP.txt > "$VIRUS"_CP.txt
    grep -F -v "$VIRUS" NeighborSumValues_CP.txt > "$VIRUS"_testCP.txt
done
