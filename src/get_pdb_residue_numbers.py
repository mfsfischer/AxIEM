#!/usr/bin/env python3

#author: Marion Fischer
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBExceptions
from Bio.SeqUtils import seq1
import sys

parser = PDBParser()

def usage():
    print("Prints PDB residue numbers - for data table curation only")
    print("python get_fasta_from_pdb_residue_numbers.py <pdb>")

    warnings.simplefilter('ignore',PDBExceptions.PDBConstructionWarning)

    if (len(sys.argv) < 1):
        usage()
        exit()

def get_residue_id_number(pdb):
    res_num = []
    residues = [r for r in structure.get_residues() if r.get_id()[0] == " "]
    for resi in residues:
         res_num.append(resi.get_id()[1])
    return res_num

structure = parser.get_structure("apdb", sys.argv[1])
res_num = get_residue_id_number(structure)
for r in res_num:
    print(r)

     
