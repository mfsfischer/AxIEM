#!/bin/bash

for PDB in $(cat pdb.list)

do
    python ../src/insert_bfactor.py -p ../pdb_structures/"$PDB".pdb -b "$PDB"_b.txt -o "$PDB"_predictions.pdb
done
