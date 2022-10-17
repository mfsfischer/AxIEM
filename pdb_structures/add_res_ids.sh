#!/bin/bash

touch res-ids.txt
for PDB in $(cat order.list)
do
    cat "$PDB"_res-ids.txt >> res-ids.txt 
done
