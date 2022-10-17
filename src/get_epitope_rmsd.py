#!/usr/bin/env python3

import numpy
import pandas
import sys
import re
from argparse import ArgumentParser
from collections import defaultdict
from io import StringIO

parser = ArgumentParser(
    """
    Returns the RMSD of two common epitopes and each epitope's centroid as an xyz coordinate and residue closest to xyz coordinates
    """
)
parser.add_argument(
    "--pdb1",
    required=True,
    dest="pdb1",
    help="""PDB of virus 1"""
)
parser.add_argument(
    "--pdb2",
    required=True,
    dest="pdb2",
    help="""PDB of virus 2"""
)
parser.add_argument(
    "--res1",
    required=True,
    dest="res1",
    help="""Input file containing two columns '<chain ID> <residue#>' for epitope of virus 1"""
)
parser.add_argument(
    "--res2",
    required=True,
    dest="res2",
    help="""Input file containing two columns '<chain ID> <residue#>' for epitope of virus 2"""
)
parser.add_argument(
    "--rmsd",
    required=False,
    default=0,
    dest="compare",
    help="""Compare RMSD of two epitopes: 1 (yes) or 0 (no); default is 0"""
)

args = parser.parse_args()

def getEpitopeCoordinates(pdb, chains, resnum):
    """
    Extracts xyz coordinates of residues listed as epitopes
    @params pdb : raw pdb file
    @params chains : chain ID of epitope residues of size (n)
    @params resnum : residue number of epitope residues of size (n)
    @returns res_coords : returns array of xyz coordinates for all epitopes residues of size (n,3)
    """
    x_coords = []
    y_coords = []
    z_coords = []
    i = 0
    with open(pdb, 'r') as f:
        pdb = f.readlines()
        for line in pdb:
            chain = line[21:22]
            number = line[22:26]
            if i < len(resnum) and chain==chains[i] and str(number)==str(resnum[i]):
                x_coords.append(float(line[30:38]))
                y_coords.append(float(line[38:46]))
                z_coords.append(float(line[46:54]))
                i += 1
                
    res_coords = [x_coords, y_coords, z_coords]
    res_coords = numpy.asarray(res_coords)
    return res_coords

def getEpitopeCenter(pdb, res_coords, n):
    """
    @params pdb : raw pdb file
    @params res_coords : xyz coords of epitope residues
    @returns center : xyz coordinate of epitope center
    @returns centralresidue : residue closest to center
    """
    center = numpy.sum(res_coords, axis=1) / n

    # Find residue nearest epitope center
    distance=1000
    centralresidue = []
    with open(pdb, 'r') as f:
        pdb = f.readlines()
        for line in pdb:
            if line[12:16]==" CA ":
                ca_distance = numpy.sqrt( (float(line[30:38])-center[0])**2 + (float(line[38:46])-center[1])**2 + (float(line[46:54])-center[2])**2 )
                if ca_distance < distance:
                    distance = ca_distance
                    centralresidue = line[21:26]
    return center, centralresidue

def getEpitopeRMSD(res_coords1, res_coords2):
    """
    @params res_coords1 : x coordinates of residues in epitope 1
    @params res_coords2 : x coordinates of residues in epitope 2
    @returns rmsd : root mean square difference of pairwise distances b/w 1 and 2
    """
    rmsd = 0
    if res_coords1.shape != res_coords2.shape:
        print("Epitope 1 has", len(res_coords1), "residues while epitope 2 has", len(res_coords2), "residues")
        print("Cannot compute rmsd- both must have same # of residues. Quitting.")
        sys.exit(-1)
    else:
        #differences = numpy.subtract( res_coords1, res_coords2 )
        #d2 = 0
        #for r in enumerate(res_coords1):
        #    d2 += (differences[r[0]][0]**2 + differences[r[0]][1]**2 + differences[r[0]][2])
        #rmsd = numpy.sqrt(d2/len(res_coords1))
        rmsd_res = 0
        for i in range(res_coords1.shape[1]):
            for j in range(res_coords1.shape[1]):
                d1 = (res_coords1[0][j] - res_coords1[0][i])**2 + (res_coords1[1][j] - res_coords1[1][i])**2 + (res_coords1[2][j] - res_coords1[2][i])**2 
                d2 = (res_coords2[0][j] - res_coords2[0][i])**2 + (res_coords2[1][j] - res_coords2[1][i])**2 + (res_coords2[2][j] - res_coords2[2][i])**2
                rmsd_res += d1-d2
            rmsd = rmsd + numpy.abs(rmsd_res/res_coords1.shape[1])
        rmsd = rmsd/res_coords1.shape[1]
        rmsd = numpy.sqrt(rmsd)
            
    return rmsd
            
            

def main():

    residues1 = pandas.read_csv(args.res1, header=None, delimiter=' ')
    chains1 = residues1.iloc[:,0].to_list()
    resnum1 = residues1.iloc[:,1].to_numpy()
    residues2 = pandas.read_csv(args.res2, header=None, delimiter=' ')
    chains2 = residues2.iloc[:,0].to_list()
    resnum2 = residues2.iloc[:,1].to_numpy()

    res_coords1 = getEpitopeCoordinates(args.pdb1, chains1, resnum1)
    res_coords2 = getEpitopeCoordinates(args.pdb2, chains2, resnum2)
    #print(res_coords1)

    center1, centralresidue1 = getEpitopeCenter(args.pdb1, res_coords1, len(resnum1))
    center2, centralresidue2 = getEpitopeCenter(args.pdb2, res_coords2, len(resnum2))
    if int(args.compare)==1:
        rmsd = getEpitopeRMSD(res_coords1, res_coords2)

    print("Protein\t\t\t\t center X\t\t center Y\t\t center Z\t\t Residue ID")
    print(args.pdb1,"\t",center1[0],"\t",center1[1],"\t",center1[2],"\t",centralresidue1)
    print(args.pdb2,"\t",center2[0],"\t",center2[1],"\t",center2[2],"\t",centralresidue2)
    print("\n")
    if int(args.compare)==1:
        print("RMSD :",rmsd)

main()
    
            
                             
            
                
            



