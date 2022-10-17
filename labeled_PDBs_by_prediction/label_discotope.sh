#!/bin/bash

for PDB in $(cat pdb.list)

do
    python ../../../../scripts/insert_bfactor.py -p ../"$PDB"_designed.pdb -b "$PDB"_discotope-b.txt -o "$PDB"_discotope.pdb
done
