#!/bin/bash

for PDB in $(cat mers.list)

do
    python ~/projects/AxIEM/src/get_pdb_residue_numbers.py "$PDB" > "$PDB"_res-ids.txt 
done
