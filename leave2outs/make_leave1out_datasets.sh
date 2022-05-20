#!/bin/bash

for VIRUS in $(cat ../virus.list)

do
    cat ../training_datasets/"$VIRUS"_test.csv | cut -d' ' -f1,4-5 > "$VIRUS"_test_REU.txt
    cat ../testing_datasets/"$VIRUS"_epi.csv | cut -d' ' -f1,4-5 > "$VIRUS"_REU.txt 
    cat ../training_datasets/"$VIRUS"_test.csv | cut -d' ' -f2,4-5 > "$VIRUS"_test_CP.txt
    cat ../testing_datasets/"$VIRUS"_epi.csv | cut -d' ' -f2,4-5 > "$VIRUS"_CP.txt 
    cat ../training_datasets/"$VIRUS"_test.csv | cut -d' ' -f3-5 > "$VIRUS"_test_NV.txt
    cat ../testing_datasets/"$VIRUS"_epi.csv | cut -d' ' -f3-5 > "$VIRUS"_NV.txt 
done
