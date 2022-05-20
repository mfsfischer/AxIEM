#!/usr/bin/env python3

#author: Marion Fischer
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBExceptions
from Bio.SeqUtils import seq1
import sys

parser = PDBParser()

def usage():
    print("Creates fasta file of amino acids from a PDB")
    print("python get_fasta_from_pdb_by_chain.py <pdb>")

    warnings.simplefilter('ignore',PDBExceptions.PDBConstructionWarning)

    if (len(sys.argv) < 1):
        usage()
        exit()

def get_three_letter(pdb):
    fasta_three = []
    residues = [r for r in structure.get_residues() if r.get_id()[0] == " "]
    for resi in residues:
        fasta_three.append(resi.get_resname())
    return fasta_three

def get_one_letter(list_of_three):
    fasta_one=[]
    for x in list_of_three:
        fasta_one.append(seq1(x))
    return fasta_one

structure = parser.get_structure("apdb", sys.argv[1])
three_letter = get_three_letter(structure)
one_letter = get_one_letter(three_letter)

print(">", sys.argv[1])
print("".join(one_letter))

