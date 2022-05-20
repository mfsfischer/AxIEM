#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    cat "$VIRUS"_test.csv | cut -d' ' -f2-5 > "$VIRUS"_test_noREU.txt
    cat "$VIRUS"_epi.csv | cut -d' ' -f2-5 > "$VIRUS"_noREU.txt 
    cat "$VIRUS"_test.csv | cut -d' ' -f1,3-5 > "$VIRUS"_test_noCP.txt
    cat "$VIRUS"_epi.csv | cut -d' ' -f1,3-5 > "$VIRUS"_noCP.txt 
    cat "$VIRUS"_test.csv | cut -d' ' -f1-2,4-5 > "$VIRUS"_test_noNV.txt
    cat "$VIRUS"_epi.csv | cut -d' ' -f1-2,4-5 > "$VIRUS"_noNV.txt 
done
