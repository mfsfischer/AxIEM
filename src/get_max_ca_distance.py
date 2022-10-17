#!/usr/bin/env python3
"""
Finds maximum C-alpha to C-alpha distance within a single protein given one or more PDB models.
More than one protein may be analyzed per run.

@Returns file containing 2 columns: <Protein Name> <Max C-alpha distance> for each unique protein name (note: labels are case sensitive)
"""

import math
import numpy
import re
from argparse import ArgumentParser
from Bio.PDB.PDBParser import PDBParser

parser = ArgumentParser()
parser.add_argument(
    "-l",
    "--pdb_list",
    dest="pdbs",
    required=True,
    help="""Input space-delimited file containing two columns: <PDB name> <PROTEIN name>. Note protein name must be identical (case-sensitive)."""
)
parser.add_argument(
    "-o",
    "--output",
    required=True,
    help="""Output space-delimited file containing two columns: <PROTEIN name> <MAX DISTANCE>."""
)

args = parser.parse_args()
pdbparser = PDBParser()

def getPDBmaxDistance(pdb):
    """
    Parse PDB (accounting for residue # >= 1000) and calculate c-alpha to c-alpha distances
    @ Returns maximum c-alpha to c-alpha distance
    """
    # Parse PDB
    pdb_parsed = []
    with open(pdb, 'r') as p:
        pdb = p.readlines()
        for line in pdb:
            # Account for numerous instances where white space no longer exists between columns
            # Residue numbers >=1000
            num_index = 22
            num_line = line[:num_index] + ' ' + line[num_index:]
            # Y coordinate (add 1 index from initial index 22)
            y_index = 39
            y_line = num_line[:y_index] + ' ' + num_line[y_index:]
            # Z coordinate (add 2 index from initial index 46)
            z_index = 48
            z_line = y_line[:z_index] + ' ' + y_line[z_index:]
            pdb_parsed.append(re.split('\s+', z_line))
    # Calculate c-alpha distances and identify maximum distance
    distance = 0.0
    max_distance = 0.0
    for first in enumerate(pdb_parsed):
        if first[1][0]=="ATOM":
            ca1=numpy.zeros(3)
            if first[1][2]=="CA":
                ca1[0]=float(first[1][6])
                ca1[1]=float(first[1][7])
                ca1[2]=float(first[1][8])
                for second in enumerate(pdb_parsed):
                    ca2=numpy.zeros(3)
                    if second[1][0]=="ATOM":
                        if second[1][2]=="CA":
                            # When chain ID is not the same
                            if second[1][4]!=first[1][4]:
                                ca2[0]=float(second[1][6])
                                ca2[1]=float(second[1][7])
                                ca2[2]=float(second[1][8])
                            # When chain ID matches, only compare residue#s > first's residue# to compute half of contact matrix
                            elif second[1][4]==first[1][4] and int(second[1][5])>int(first[1][5]):
                                ca2[0]=float(second[1][6])
                                ca2[1]=float(second[1][7])
                                ca2[2]=float(second[1][8])
                    distance = numpy.sqrt( (ca2[0]-ca1[0])**2 + (ca2[1]-ca1[1])**2 + (ca2[2]-ca1[1])**2 )
                    if distance > max_distance:
                        max_distance = distance
    return max_distance

def main():
    """
    Stores max distance for each unique protein name as dictionary containing max c-alpha to c-alpha distances.
    Goes through each PDB file and either creates new key value for a unique protein name or checks to find 
    maximum c-alpha distance > than previous PDB model(s)
    """

    # Read in pdbs to parsed
    names = []
    with open(args.pdbs, 'r') as p:
        lines = p.readlines()
        for line in lines:
            names.append(line.split())

    # Iterate through each PDB to calculate c-alpha to c-alpha distance values to find max distance for each pdb
    # Either update or create dict key and value of protein name and max distance value
    protein_distances = {}
    for pdb in enumerate(names):
        print("Get max ca-ca distances within",pdb[1][0],"for protein",pdb[1][1])
        # Avoid using Biopython because it is super slow in parsing PDB files for contact matrices
        pdb_max_distance = getPDBmaxDistance(pdb[1][0])
        if pdb[1][1] in protein_distances:
            if protein_distances[pdb[1][1]] < pdb_max_distance:
                protein_distances[pdb[1][1]] = pdb_max_distance
        else:
            protein_distances[pdb[1][1]] = pdb_max_distance
                
    with open(args.output, 'w') as o:
        for protein, distance in protein_distances.items():
            o.write(protein + " " + str(distance) + "\n")
        
main()
